#ifndef __CAMERA_H__
#define __CAMERA_H__

#include <stdbool.h>

void camera_enable_dma (bool enable);
void camera_set_front_buffer (int *buff);
void camera_set_back_buffer (int *buff);
void camera_swap_buffers(void);

#endif /* __CAMERA_H__ */
