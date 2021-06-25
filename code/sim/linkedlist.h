#pragma once

#include <stdlib.h>
typedef struct cons Cons;

struct cons {
    Cons* next;
    void* value;
};

Cons* cons_new(void* v);
void cons_push_inc(Cons** c, void* v);
int cons_size(Cons* c);
Cons* cons_end(Cons* c);
void cons_free_all(Cons* c);
void cons_free_list(Cons* c);
