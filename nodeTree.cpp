//
// Created by brianpage on 5/29/17.
//

#include "nodeTree.h"

// Height of AVL Tree
int nodeTree::height(avl_node *temp){
    int h = 0;
    if (temp != NULL)
    {
        int l_height = height (temp->left);
        int r_height = height (temp->right);
        int max_height = max (l_height, r_height);
        h = max_height + 1;
    }
    return h;
}

// Height Difference
int nodeTree::diff(avl_node *temp){
    int l_height = height (temp->left);
    int r_height = height (temp->right);
    int b_factor= l_height - r_height;
    return b_factor;
}

// Right- Right Rotation
avl_node *nodeTree::rr_rotation(avl_node *parent){
    avl_node *temp;
    temp = parent->right;
    parent->right = temp->left;
    temp->left = parent;
    return temp;
}

// Left- Left Rotation
avl_node *nodeTree::ll_rotation(avl_node *parent){
    avl_node *temp;
    temp = parent->left;
    parent->left = temp->right;
    temp->right = parent;
    return temp;
}

// Left - Right Rotation
avl_node *nodeTree::lr_rotation(avl_node *parent){
    avl_node *temp;
    temp = parent->left;
    parent->left = rr_rotation (temp);
    return ll_rotation (parent);
}

// Right- Left Rotation
avl_node *nodeTree::rl_rotation(avl_node *parent){
    avl_node *temp;
    temp = parent->right;
    parent->right = ll_rotation (temp);
    return rr_rotation (parent);
}

// Balancing AVL Tree
avl_node *nodeTree::balance(avl_node *temp){
    int bal_factor = diff (temp);
    if (bal_factor > 1)
    {
        if (diff (temp->left) > 0)
            temp = ll_rotation (temp);
        else
            temp = lr_rotation (temp);
    }
    else if (bal_factor < -1)
    {
        if (diff (temp->right) > 0)
            temp = rl_rotation (temp);
        else
            temp = rr_rotation (temp);
    }
    return temp;
}

// Insert Element into the tree
avl_node *nodeTree::insert(int value, int processNumber){
    if (root == NULL){
        root = new avl_node;
        root->assignedElementCount = value;
        root->processId = processNumber;
        root->left = NULL;
        root->right = NULL;
        return root;
    } else {
        avl_node *temp = root;
        while (temp->left != NULL){
            temp = temp->left;
        }

        temp = new avl_node;
        temp->assignedElementCount = value;
        temp->processId = processNumber;
        temp->left = NULL;
        temp->right = NULL;

        root = balance (temp);

        return temp;
    }
}

// Find node with lowest assigned element count
avl_node *nodeTree::findLowest(avl_node *tree){
    //std::cout << "inside findLowest" << std::endl;
    avl_node *temp = tree;
    //std::cout << "findLowest-2" << std::endl;
    while (temp->left != NULL) {
        //std::cout << "findLowest-3" << std::endl;
        temp = temp->left;
    }
    //std::cout << "findLowest-4" << std::endl;
    return temp;
}

avl_node *nodeTree::assignElements(avl_node *tree, int numElementsToAssign, int &nodeAssignment){
    //std::cout << "--" << std::endl;
    avl_node * nodeToUpdate = findLowest(tree);
    //std::cout << "**" << std::endl;
    nodeToUpdate->assignedElementCount = numElementsToAssign;
    //std::cout << "@@" << std::endl;
    root = balance(tree);
    //std::cout << "%%" << std::endl;
    nodeAssignment = nodeToUpdate->processId;
    //std::cout << "&&" << std::endl;
    return root;
}

// Display AVL Tree
void nodeTree::display(avl_node *ptr, int level){
    int i;
    if (ptr!=NULL)
    {
        display(ptr->right, level + 1);
        printf("n");
        if (ptr == root)
            cout<<"Root -> ";
        for (i = 0; i < level && ptr != root; i++)
            cout<<"        ";
        cout<< ptr->assignedElementCount;
        display(ptr->left, level + 1);
    }
}

// Inorder Traversal of AVL Tree
void nodeTree::inorder(avl_node *tree){
    if (tree == NULL)
        return;
    inorder (tree->left);
    cout<< tree->assignedElementCount <<"  ";
    inorder (tree->right);
}

// Preorder Traversal of AVL Tree
void nodeTree::preorder(avl_node *tree){
    if (tree == NULL)
        return;
    cout<< tree->assignedElementCount <<"  ";
    preorder (tree->left);
    preorder (tree->right);

}

//  Postorder Traversal of AVL Tree
void nodeTree::postorder(avl_node *tree){
    if (tree == NULL)
        return;
    postorder ( tree ->left );
    postorder ( tree ->right );
    cout<< tree->assignedElementCount <<"  ";
}

void nodeTree::freeTree(){
    avl_node *temp = root;

    //
    // Make destructor delete all nodes from tree!
    //
}