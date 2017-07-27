//
// Created by brianpage on 5/29/17.
//

#ifndef DISTRUBUTED_SPMV_NODETREE_H
#define DISTRUBUTED_SPMV_NODETREE_H

#include<iostream>
#include<cstdio>
#include<sstream>
#include<algorithm>

using namespace std;

struct avl_node {
    int assignedElementCount;
    int processId;
    struct avl_node *left;
    struct avl_node *right;
};

class nodeTree {
    public:
        int height(avl_node *);
        int diff(avl_node *);
        avl_node *root = NULL;

        avl_node *rr_rotation(avl_node *);
        avl_node *ll_rotation(avl_node *);
        avl_node *lr_rotation(avl_node *);
        avl_node *rl_rotation(avl_node *);
        avl_node *balance(avl_node *);
        avl_node *insert(int value, int processNumber);
        avl_node *findLowest(avl_node *);
        avl_node *assignElements(avl_node *root, int numElementsToAssign, int &nodeAssignment);

        void display(avl_node *, int);
        void inorder(avl_node *);
        void preorder(avl_node *);
        void postorder(avl_node *);
        void freeTree();

        nodeTree() {
            root = NULL;
            //root->processId = 0;
        }

        ~nodeTree(){
            freeTree();
        }
};

#endif //DISTRUBUTED_SPMV_NODETREE_H
