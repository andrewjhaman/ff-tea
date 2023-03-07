
#scope_export

#import "Basic";

#scope_export 
image : Image(u8);
Image :: struct ($T : Type) {
    width  : s32;
    height : s32;
    channel_count : s32 = 1;
    data : [] T;
}


fft :: inline (values : [] $T) -> [] Complex  {
    //values.count must be a power of 2.
    return __fft_impl(values);
}

ifft :: inline (values : [] $T) -> [] Complex {
    //values.count must be a power of 2
    result := __fft_impl(values, is_inverse = true);
    recip := 1.0 / cast(float)result.count;
    for *v : result {
        <<v *= recip;
    }
    return result;
}


image_fft :: (image : Image) -> Image(Complex) {
    input := to_grayscale(image, float);
    result : Image(Complex);
    result.data = NewArray(image.width * image.height, Complex);
    result.width = image.width;
    result.height = image.height;
    result.channel_count = 1;

    intermediate := NewArray(image.width * image.height, Complex, allocator=temp);
    for x: 0..image.width-1 {
        column := NewArray(image.height, float, allocator=temp);
        for y : 0..image.height-1 {
            column[y] = input.data[x + y * image.width];
        }
        column_fft := fft(column);
        for y : 0..image.height-1 {
            intermediate[x + y * image.width] = column_fft[y];
        }
    }

    for y : 0..image.height-1 {
        row : [] Complex;
        row.data = intermediate.data + y * image.width;
        row.count = image.width;
        row_fft := fft(row);

        memcpy(result.data.data + (y * image.width), row_fft.data, row_fft.count * size_of(type_of(row_fft[0])));
    }

    return result;
}

image_ifft :: (image : Image(Complex)) -> Image(float) {
    input := image;
    result := NewImage(image.width, image.height);
    
    intermediate := NewArray(image.width * image.height, Complex, allocator=temp);
    for x: 0..image.width-1 {
        column := NewArray(image.height, Complex, allocator=temp);
        for y : 0..image.height-1 {
            column[y] = input.data[x + y * image.width];
        }
        column_fft := ifft(column);
        for y : 0..image.height-1 {
            intermediate[x + y * image.width] = column_fft[y];
        }
    }

    for y : 0..image.height-1 {
        row : [] Complex;
        row.data = intermediate.data + y * image.width;
        row.count = image.width;
        row_fft := ifft(row);

        for v : row_fft {
            result.data[it_index + (y * image.width)] = v.x;
        }
    }

    return result;
}




NewImage :: (width : s32, height : s32, channel_count : s32 = 1, $data_type : Type = float) -> Image(data_type) {
    image : Image(data_type);
    image.width = width;
    image.height = height;
    image.channel_count = channel_count;
    image.data = NewArray(width * height, data_type);
    return image;
}

NewImage :: (image : Image, channel_count : s32 = 1, $data_type : Type) -> Image(data_type) {
    return NewImage(image.width, image.height, channel_count, data_type);
}

NewImage :: (image : Image($T)) -> Image(T) {
    return NewImage(image.width, image.height, image.channel_count, T);
}

copy_image :: (image : Image($T)) -> Image(T) {
    result := image;
    result.data = array_copy(image.data);
    return result;
}

convert_type :: (image : Image, $type : Type, scale_factor := 1.0) -> Image(type) {
    result : Image(type);
    result.width = image.width;
    result.height = image.height;
    result.channel_count = image.channel_count;
    result.data = NewArray(image.width * image.height, type);

    for v : image.data {
        result.data[it_index] = cast(type)(v * scale_factor);
    }
    return result;
}

to_grayscale :: (image : Image, $type : Type = float) -> Image(type) {
    #if type_of(image.data[0]) == type
        if image.channel_count == 1 then return copy_image(image);

    result : Image(type);
    result.width = image.width;
    result.height = image.height;
    result.channel_count = 1;
    result.data = NewArray(result.width * result.height, type);

    for y : 0..result.height-1 {
        for x : 0..result.width - 1 {
            index := x + y*image.width;
            index *= image.channel_count;

            //Just do a dumb average right now.
            sum := 0.0;
            for 0..image.channel_count-1 { sum += cast(float)image.data[index+it];}
            sum /= cast(float)image.channel_count;
            result.data[x + y*result.width] = cast(type)sum;
        }
    }
    return result;
}


get_real_part :: (image : Image(Complex)) -> Image(float) {
    result : Image(float);
    result.data = NewArray(image.width * image.height, float);
    result.width = image.width;
    result.height = image.height;
    result.channel_count = 1;

    for v : image.data {
        result.data[it_index] = v.x;
    }
    return result;
}

get_imag_part :: (image : Image(Complex)) -> Image(float) {
    result : Image(float);
    result.data = NewArray(image.width * image.height, float);
    result.width = image.width;
    result.height = image.height;
    result.channel_count = 1;

    for v : image.data {
        result.data[it_index] = v.y;
    }
    return result;
}

get_magnitude :: (image : Image(Complex)) -> Image(float) {
    result : Image(float);
    result.data = NewArray(image.width * image.height, float);
    result.width = image.width;
    result.height = image.height;
    result.channel_count = 1;

    for v : image.data {
        result.data[it_index] = length(v);
    }
    return result;
}


shift_topleft_to_centre :: (image : Image($T)) -> Image(T) {
    result := NewImage(image);
    for y : 0..image.height-1 {
        for x : 0..image.width-1 {
            read_x := x + image.width / 2;
            read_y := y + image.height / 2;
            read_x %= image.width;
            read_y %= image.height;
            result.data[x + y*image.width] = image.data[read_x + read_y * image.width];
        }
    }
    return result;
}

to_fft_canonical_form :: (image : Image(Complex)) -> Image(float) {
    magnitude := get_magnitude(image);

    maxitude := -FLOAT32_INFINITY;
    for v : magnitude.data {
        maxitude = max(v, maxitude);
    }

    c := 1.0 / log(1 + maxitude);

    for *v : magnitude.data {
        <<v = c * log(1.0 + <<v);
    }

    shifted := shift_topleft_to_centre(magnitude);

    return shifted;
}


make_complex :: (real : float, imaginary : float) -> Complex {
    result : Complex = ---;
    result.x = real;
    result.y = imaginary;
    return result;
}

exp :: (c : Complex) -> Complex {
    result : Complex;

    //e^(a+bi)=(e^a)*e^(bi)=(e^a)(cos(b) + i*sin(b))

    e_a := exp(c.x);
    e_b := make_complex(cos(c.y), sin(c.y));

    result = e_a * e_b;

    return result;
}


slow_dft :: (values : [] $T) -> [] Complex {
    //Works for now, but doesn't support inverses yet.
    result := NewArray(values.count, Complex);

    i := make_complex(0, 1);
    for k : 0..values.count-1 {
        for v : values {
            c := exp(make_complex(0, -TAU / cast(float)values.count * cast(float)k * cast(float)it_index));
            t := v * c;
            result[k] += t;
        }
    }

    return result;
}




#scope_file
#import "Math";


example_usage :: () {
    values : [1<<3]float;

    period : float = 128.0;

    for *v : values {
        <<v = cos(it_index * TAU / cast(float)values.count);
        <<v += cos(it_index * TAU / cast(float)values.count * 31);
    }   

    fft_values := fft(values);
    reconstructed := ifft(fft_values);
}




__fft_impl :: (values : [] $T, $is_inverse := false) -> [] Complex {
    N := values.count;

    if N <= 1 {//TODO(Andrew): Determine the best value here
        result := NewArray(1, Complex, allocator=temp);
        #if T == Complex then
            result[0] = values[0];
        else
            result[0] = .{values[0], 0.0};
        return result;
        //return slow_dft(values);
    }

    result := NewArray(values.count, Complex);

    even := NewArray(N / 2, T, allocator=temp);
    odd := NewArray(N / 2, T, allocator=temp);

    for 0..values.count/2-1 {
        even[it] = values[it*2  ];
        odd [it] = values[it*2+1];
    }

    even_c := __fft_impl(even, is_inverse);
    odd_c := __fft_impl(odd, is_inverse);

    for 0..(N/2)-1 {
        //Can probably speed up this exp by using fact that real part is 0
        inverse_factor := ifx is_inverse then -1.0 else 1.0;
        factor : Complex = exp(make_complex(0, -TAU * cast(float)it / cast(float)N * inverse_factor));
        odd_term : Complex = factor * odd_c[it];

        result[it] = even_c[it] + odd_term;
        result[it+N/2] = even_c[it] - odd_term;
    }

    return result;
}