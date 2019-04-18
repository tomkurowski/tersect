/*  heap.c

    Copyright (C) 2019 Cranfield University

    Author: Tomasz Kurowski <t.j.kurowski@cranfield.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE. */

#include "heap.h"

#include <stdlib.h>

// Get index of left/right child elements for a heap element
#define LEFT_CHILD(index) (index << 1) + 1
#define RIGHT_CHILD(index) (index << 1) + 2
#define PARENT(index) (index - 1) >> 1

/* Example comparator (for ints) */
/* int intcomp(const void *a, const void *b)
{
    return *(const int *)a - *(const int *)b;
}*/

Heap *init_heap(int max_size, int(*comparator)(const void* a, const void* b))
{
    Heap *heap = malloc(sizeof *heap);
    heap->size = 0;
    heap->max_size = max_size;
    heap->array = malloc(heap->max_size * sizeof *(heap->array));
    heap->comparator = comparator;
    return heap;
}

void free_heap(Heap *heap)
{
    free(heap->array);
    free(heap);
}

static inline void _swap(void **a, void **b)
{
    int *tmp = *a;
    *a = *b;
    *b = tmp;
}

static inline void _sift_up(Heap *heap)
{
    int position = heap->size - 1;
    int parent_position = PARENT(position);
    while (position && heap->comparator(heap->array[position],
                                        heap->array[parent_position])
                                        < 0) {
        _swap(&heap->array[position], &heap->array[parent_position]);
        position = parent_position;
        parent_position = PARENT(position);
    }
}

void heap_push(Heap *heap, void *value)
{
    heap->array[heap->size++] = value;
    _sift_up(heap);
}

static inline void _sift_down(Heap *heap)
{
    int position = 0;
    int child_left;
    int child_right;
    int larger_child;
    while ((child_left = LEFT_CHILD(position)) < heap->size) { // while left child exists
        child_right = RIGHT_CHILD(position);
        if (child_right >= heap->size // right child does not exist
            || heap->comparator(heap->array[child_left],
                                heap->array[child_right])
                                < 0) {
            larger_child = child_left;
        } else {
            larger_child = child_right;
        }
        if (heap->comparator(heap->array[larger_child],
                             heap->array[position])
                             < 0) {
            _swap(&heap->array[position], &heap->array[larger_child]);
            position = larger_child;
        } else {
            break;
        }
    }
}

void sift_down(Heap *heap)
{
    _sift_down(heap);
}

void sift_up(Heap *heap)
{
    _sift_up(heap);
}

void clear_heap(Heap *heap)
{
    // No need to actually clear the array
    heap->size = 0;
}

void *heap_pop(Heap *heap)
{
    void *output = heap->array[0];
    heap->array[0] = heap->array[--(heap->size)];
    _sift_down(heap);
    return output;
}

void *heap_peek(const Heap *heap)
{
    return heap->array[0];
}
