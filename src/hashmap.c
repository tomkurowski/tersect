/*  hashmap.c

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

#include "hashmap.h"

#include <stdlib.h>
#include <string.h>

HashMap *init_hashmap(const uint32_t bucket_count)
{
    HashMap *hm = malloc(sizeof *hm);
    hm->bucket_count = bucket_count;
    hm->count = 0;
    hm->buckets = calloc(bucket_count, sizeof *(hm->buckets));
    hm->elements = calloc(bucket_count, sizeof *(hm->elements)); // this should be e.g. half of bucket count, to be doubled when load factor exceeds 0.5
    return hm;
}

static uint32_t hash_func(const char *key)
{
    uint32_t hash = 0;
    char c;
    while ((c = (char)*key++) != 0) {
        hash += c;
        hash *= 0x9e3779b1; // Knuth
    }
    return hash;
}

void hashmap_insert(HashMap *hm, const char *key, void *value)
{
    Bucket *bucket = &(hm->buckets[hash_func(key) % hm->bucket_count]);
    HashElement **element = &(bucket->first);
    for (uint32_t i = 0; i < bucket->count; ++i) {
        if (!strcmp((*element)->key, key)) {
            // Overwrite old value
            (*element)->value = value;
            return;
        }
        element = &((*element)->next);
    }
    // Adding new element at end of bucket chain (and to the list of keys)
    (*element) = &(hm->elements[hm->count]);
    // Make this some sort of string pool instead of an alloc & make sure it's freed
    (*element)->key = malloc(strlen(key) + 1);
    strcpy((*element)->key, key);
    (*element)->value = value;
    ++(hm->count);
    ++(bucket->count);
}

void *hashmap_get(const HashMap *hm, const char *key)
{
    Bucket *bucket = &(hm->buckets[hash_func(key) % hm->bucket_count]);
    HashElement *element = bucket->first;
    while (element) {
        if (!strcmp(element->key, key)) {
            return element->value;
        }
        element = element->next;
    }
    return NULL;
}

HashIterator hashmap_iterator(const HashMap *hm)
{
    HashIterator it;
    it._nextindex = 0;
    it._hm = hm;
    return it;
}

void hashmap_iterator_next(HashIterator* it)
{
    HashElement *element;
    if (it->_hm->count > it->_nextindex) {
        element = &(it->_hm->elements[it->_nextindex++]);
        it->key = element->key;
        it->value = element->value;
    } else {
        it->key = NULL;
        it->value = NULL;
    }
}

static inline void free_keys(HashMap *hm)
{
    HashIterator it = hashmap_iterator(hm);
    hashmap_iterator_next(&it);
    while (it.key) {
        free(it.key);
        hashmap_iterator_next(&it);
    }
}

void free_hashmap(HashMap *hm)
{
    free_keys(hm);
    free(hm->buckets);
    free(hm->elements);
    free(hm);
}
