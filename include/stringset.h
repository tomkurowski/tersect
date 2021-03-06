/*  stringset.h

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

#ifndef STRINGSET_H
#define STRINGSET_H

#include "hashmap.h"

#include <stdbool.h>

struct StringSet {
    HashMap *map;
};

struct StringSetIterator {
    HashIterator hm_iterator;
};

/**
 * Very basic unordered string set implementation using struct HashMap.
 */

struct StringSet *init_stringset();
void free_stringset(struct StringSet *set);

void stringset_add(const struct StringSet *set, char *str);
bool stringset_contains(const struct StringSet *set, const char *str);
uint32_t stringset_size(const struct StringSet *set);

struct StringSetIterator stringset_iterator(const struct StringSet *set);
char *stringset_iterator_next(struct StringSetIterator* it);

#endif
