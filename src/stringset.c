/*  stringset.c

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

#include "stringset.h"

#include <stdlib.h>

#define DEFAULT_STRINGSET_BUCKET_SIZE 10000

struct StringSet *init_stringset()
{
    struct StringSet *set = malloc(sizeof *set);
    set->map = init_hashmap(DEFAULT_STRINGSET_BUCKET_SIZE);
    return set;
}

void free_stringset(struct StringSet *set)
{
    free_hashmap(set->map);
    free(set);
}

void stringset_add(const struct StringSet *set, char *str)
{
    static int dummy;
    hashmap_insert(set->map, str, &dummy);
}

bool stringset_contains(const struct StringSet *set, char *str)
{
    return hashmap_get(set->map, str) != NULL;
}

uint32_t stringset_size(const struct StringSet *set)
{
    return set->map->count;
}

struct StringSetIterator stringset_iterator(const struct StringSet *set)
{
    struct StringSetIterator it;
    it.hm_iterator = hashmap_iterator(set->map);
    return it;
}

char *stringset_iterator_next(struct StringSetIterator* it)
{
    hashmap_iterator_next(&it->hm_iterator);
    return it->hm_iterator.key;
}
