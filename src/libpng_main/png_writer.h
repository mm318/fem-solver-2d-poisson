#ifndef _PNG_WRITER_H_
#define _PNG_WRITER_H_

#include <stdint.h>


/* A coloured pixel. */
typedef struct {
    uint8_t red;
    uint8_t green;
    uint8_t blue;
} pixel_t;

/* A picture. */
typedef struct {
    int height;
    int width;
    pixel_t *pixels;
    // uint8_t *pixels;
} bitmap_t;

/* structure to store PNG image bytes */
typedef struct {
    char *buffer;
    size_t size;
} FileBuffer;


bitmap_t* new_bitmap(const int height, const int width);

int set_pixel(bitmap_t * const bmp, const int x, const int y,
              const pixel_t * const color);

int fill_circle(bitmap_t * const bmp, const int x, const int y, const int radius,
                const pixel_t * const color);

int draw_line(bitmap_t * const bmp, const int x_a, const int y_a, const int x_b,
              const int y_b, const int thickness, const pixel_t * const color);

/* Write "bitmap" to a PNG file specified by "path"; returns 0 on
   success, non-zero on error. */
int write_bitmap(bitmap_t * const bitmap, const char * const path);

FileBuffer* get_png(bitmap_t * const bitmap);

int delete_bitmap(bitmap_t * const bmp);

int delete_filebuf(FileBuffer * const buf);


#endif
