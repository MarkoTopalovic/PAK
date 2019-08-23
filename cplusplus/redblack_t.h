#ifndef _RED_BLACK_H_
#define _RED_BLACK_H_

#include <stdio.h>
#include <stdint.h>
#include <malloc.h>

typedef int Tkey;
typedef double Tinfo;

typedef enum {red,black} COLOR;

typedef struct RBNode_t_tag
{
	Tkey keyRow;
	Tkey keyCol;
	Tinfo info;
	struct RBNode_t_tag *father;
	struct RBNode_t_tag *leftson;
	struct RBNode_t_tag *rightson;
	COLOR color;
} RBNode_t;

typedef void (*rb_callback_func_t)(RBNode_t * pNode);

void CountNode_callback(RBNode_t * pNode);
void BrisiNode_callback(RBNode_t * pNod);

int64_t RBCount(RBNode_t * pNode);
void RBInit();
void clean(RBNode_t * p);
void insert(Tkey kRow, Tkey kCol, Tinfo inf);
void insert2(RBNode_t * q);
RBNode_t * search(RBNode_t *p, Tkey kRow, Tkey kCol);
void inorder(RBNode_t *p);
void preorder(RBNode_t *p);
void postorder(RBNode_t *p);
RBNode_t * minimum(RBNode_t *p);
RBNode_t * maximum(RBNode_t *p);
RBNode_t * successor(RBNode_t *p);
RBNode_t * predecessor(RBNode_t *p);
void leftrotate(RBNode_t *x);
void rightrotate(RBNode_t *x);
void rb_insert(RBNode_t *x);
void rb_delete(RBNode_t *z);
void rb_delete_fixup(RBNode_t *x);

#endif


