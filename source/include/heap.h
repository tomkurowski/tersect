/*  heap.h

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

#ifndef HEAP_H
#define HEAP_H

typedef struct Heap Heap;

struct Heap {
    int size;
    int max_size;
    int(*comparator)(const void* a, const void* b);
    void **array;
};

Heap *init_heap(int max_size, int(*comparator)(const void* a, const void* b));
void heap_push(Heap *heap, void *value);
void *heap_pop(Heap *heap);
void *heap_peek(const Heap *heap);
void sift_down(Heap *heap);
void sift_up(Heap *heap);
void clear_heap(Heap *heap);
void free_heap(Heap *heap);

#endif
