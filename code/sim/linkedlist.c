#include "linkedlist.h"

Cons* cons_new(void* v) {
    Cons* c = (Cons*) malloc(sizeof(Cons));
    c->next = NULL;
    c->value = v;
    return c;
}
void cons_push_inc(Cons** c, void* v) {
    if (*c == NULL) {
        *c = cons_new(v);
    }
    else{
        (*c)->next = cons_new(v);
        *c = (*c)->next;
    }
}
int cons_size(Cons* c) {
    int size = 0;
    while (c->next != NULL) {
        size++;
        c = c->next;
    }
    return size;
}
Cons* cons_end(Cons* c) {
    while (c->next != NULL) {
        c = c->next;
    }
    return c;
}
void cons_free_all(Cons* c) {
    if (c->next != NULL) {
        cons_free_all(c->next);
    }
    free(c->value);
    free(c);
}
void cons_free_list(Cons* c) {
    if (c->next != NULL) {
        cons_free_list(c->next);
    }
    free(c);
}
