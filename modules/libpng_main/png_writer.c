#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "./libpng/png.h"
#include "png_writer.h"


bitmap_t* new_bitmap(const int height, const int width)
{
    bitmap_t *bmp;
    bmp = (bitmap_t*) malloc(sizeof(bitmap_t));
    bmp->pixels = (pixel_t*) malloc(height*width*sizeof(pixel_t));
    bmp->height = height;
    bmp->width = width;
    return bmp;
}


/* Given "bitmap", this returns the pixel of bitmap at the point ("x", "y"). */
static pixel_t * pixel_at (bitmap_t * bitmap, int x, int y)
{
    return bitmap->pixels + bitmap->width * y + x;
}


int set_pixel(bitmap_t * const bmp, const int x, const int y,
              const pixel_t * const color)
{
    if(y < 0 || y >= bmp->height || x < 0 || x >= bmp->width) {
        return -1;
    }

    int pixel_index;

    pixel_index = y*bmp->width+x;
    bmp->pixels[pixel_index].red = color->red;
    bmp->pixels[pixel_index].green = color->green;
    bmp->pixels[pixel_index].blue = color->blue;

    return 0;
}

int fill_circle(bitmap_t * const bmp, const int x, const int y, const int radius,
                const pixel_t * const color)
{
    // actual radius includes the pixel (x, y) itself, so take int radius - 1
    int i, j;
    int y_offset;

    // int r = radius-1;
    for(i=-(radius-1); i<radius; i++) {
        y_offset = (int) (sqrt((double)(radius*radius - i*i)) + 0.5);
        for(j=-(y_offset-1); j<y_offset; j++) {
            set_pixel(bmp, x+i, y+j, color);
        }
    }

    return 0;
}

int draw_line(bitmap_t * const bmp, const int x_a, const int y_a, const int x_b,
              const int y_b, const int thickness, const pixel_t * const color)
{
    int i, j;
    int x, y;
    int dx, dy;

    dx = x_b - x_a;
    dy = y_b - y_a;

    i = 0;
    if(abs(dy) > abs(dx)) {
        while(1) {
            x = x_a + (int)(dx*((double)i)/((double)dy) + 0.5);
            y = y_a + i;
            set_pixel(bmp, x, y, color);
            // for(j=1;j<thickness;j++) {
            // set_pixel(bmp, x, y-j, color);
            // set_pixel(bmp, x, y+j, color);
            // }
            for(j=1; j<thickness; j++) {
                set_pixel(bmp, x-j, y, color);
                set_pixel(bmp, x+j, y, color);
            }

            if(i < dy) {
                i++;
            } else if (i > dy) {
                i--;
            } else {
                break;
            }
        }
    } else {
        while(1) {
            x = x_a + i;
            y = y_a + (int)(dy*((double)i)/((double)dx) + 0.5);
            set_pixel(bmp, x, y, color);
            for(j=1; j<thickness; j++) {
                set_pixel(bmp, x, y-j, color);
                set_pixel(bmp, x, y+j, color);
            }
            // for(j=1;j<thickness;j++) {
            // set_pixel(bmp, x-j, y, color);
            // set_pixel(bmp, x+j, y, color);
            // }

            if(i < dx) {
                i++;
            } else if (i > dx) {
                i--;
            } else {
                break;
            }
        }
    }

    fill_circle(bmp, x_a, y_a, thickness, color);
    fill_circle(bmp, x_b, y_b, thickness, color);

    return 0;
}

void my_png_write_data(png_structp png_ptr, png_bytep data, png_size_t length)
{
    FileBuffer *p = (FileBuffer*) png_get_io_ptr(png_ptr);
    size_t nsize = p->size + length;

    /* allocate or grow buffer */
    p->buffer = realloc(p->buffer, nsize);
    if(!p->buffer)
        png_error(png_ptr, "Write Error");

    /* copy new bytes to end of buffer */
    memcpy(p->buffer + p->size, data, length);
    p->size += length;
}

/* This is optional but included to show how png_set_write_fn() is called */
// void my_png_flush(png_structp png_ptr)
// {
// }

FileBuffer* get_png(bitmap_t * const bitmap)
{
    size_t x, y;
    png_structp png_ptr = NULL;
    png_infop info_ptr = NULL;
    png_byte **row_pointers = NULL;

    /* The following number is set by trial and error only. I cannot
       see where it it is documented in the libpng manual.
    */
    static const int pixel_size = 3;
    static const int depth = 8;

    FileBuffer *buf = (FileBuffer*) malloc(sizeof(FileBuffer));
    buf->buffer = NULL;
    buf->size = 0;

    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (png_ptr == NULL) {
        goto png_create_write_struct_failed;
    }

    info_ptr = png_create_info_struct(png_ptr);
    if (info_ptr == NULL) {
        goto png_create_info_struct_failed;
    }

    /* Set up error handling. */
    if (setjmp(png_jmpbuf(png_ptr))) {
        goto png_failure;
    }

    /* Set image attributes. */
    png_set_IHDR (png_ptr,
                  info_ptr,
                  bitmap->width,
                  bitmap->height,
                  depth,
                  PNG_COLOR_TYPE_RGB,
                  PNG_INTERLACE_NONE,
                  PNG_COMPRESSION_TYPE_DEFAULT,
                  PNG_FILTER_TYPE_DEFAULT);

    /* Initialize rows of PNG. */
    row_pointers = png_malloc (png_ptr, bitmap->height * sizeof (png_byte *));
    for (y = 0; y < bitmap->height; ++y) {
        png_byte *row =
            png_malloc (png_ptr, sizeof(uint8_t) * bitmap->width * pixel_size);
        row_pointers[y] = row;
        for (x = 0; x < bitmap->width; ++x) {
            pixel_t *pixel = pixel_at(bitmap, x, y);
            *row++ = pixel->red;
            *row++ = pixel->green;
            *row++ = pixel->blue;
        }
    }

    png_set_rows(png_ptr, info_ptr, row_pointers);

    /* if my_png_flush() is not needed, change the arg to NULL */
    png_set_write_fn(png_ptr, buf, my_png_write_data, NULL);
    // png_set_write_fn(png_ptr, buf, my_png_write_data, my_png_flush);

    png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);

    for (y = 0; y < bitmap->height; y++) {
        png_free(png_ptr, row_pointers[y]);
    }
    png_free(png_ptr, row_pointers);

png_failure:
png_create_info_struct_failed:
    png_destroy_write_struct (&png_ptr, &info_ptr);
png_create_write_struct_failed:
    return buf;
}


/* Write "bitmap" to a PNG file specified by "path"; returns 0 on
   success, non-zero on error. */
int write_bitmap(bitmap_t * const bitmap, const char * const path)
{
    FileBuffer *buf;
    FILE *fp;

    buf = get_png(bitmap);
    fp = fopen (path, "wb");
    if (!fp) {
        fprintf(stderr,"failed to open \"%s\" for write!\n", path);
        return -1;
    }

    fwrite(buf->buffer, sizeof(char), buf->size, fp);
    fclose(fp);
    delete_filebuf(buf);

    return 0;
}


int delete_bitmap(bitmap_t * const bmp)
{
    free(bmp->pixels);
    free(bmp);
    return 0;
}

int delete_filebuf(FileBuffer * const buf)
{
    free(buf->buffer);
    free(buf);
    return 0;
}
