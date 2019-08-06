/*  hashmap.h

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

#ifndef HASHMAP_H
#define HASHMAP_H

#include <stdint.h>

// Forward declarations
typedef struct HashMap HashMap;
typedef struct Bucket Bucket;
typedef struct HashElement HashElement;
typedef struct HashIterator HashIterator;

struct HashMap {
    uint32_t count;
    uint32_t bucket_count;
    HashElement *elements;
    Bucket *buckets;
};

struct Bucket {
    uint32_t count;
    HashElement *first;
};

struct HashElement {
    char *key;
    void *value;
    HashElement *next;
};

struct HashIterator {
    // Public
    char *key;
    void *value;
    // Private
    const HashMap *_hm;
    uint32_t _nextindex;
};

HashIterator hashmap_iterator(const HashMap *hm);
void hashmap_iterator_next(HashIterator* it);

HashMap *init_hashmap(uint32_t bucket_count);
void hashmap_insert(HashMap *hm, const char *key, void *value);
void *hashmap_get(const HashMap *hm, const char *key);
void free_hashmap(HashMap *hm);

#endif
