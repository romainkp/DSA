
#include <R.h> 
#include "dsa_tree.h"

#define __DEBUG_TREE__ 0

p_node_t* make_new_node(unsigned long name, double rss) 
{
    /** need to allocate room for a new p_node, and correspondingly
	make some room for its children. **/
    p_node_t* root = (p_node_t*) Calloc(1, p_node_t);
    root->children = (p_node_t*) Calloc(BASE_SIZE, p_node_t);
    root->current_children = 0;
    root->max_children = BASE_SIZE;
    root->name = name;
    root->rss = rss; 

    return root;
}



/**
 * recurse while we are equal if we have reached 'len' then the 
 * node is in the tree if the RSS is filled in then we have fit the 
 * model, otherwise we haven't fit the model.
 */
p_node_t* get_node (unsigned long* probe, int len, p_node_t* tree) {
    int i = 0, j = 0, inner_match = 0, outer_match = 0;
    p_node_t* current_node = tree;

    if (__DEBUG_TREE__) {
	Rprintf("in get_node (len = %d): ", len);
	for (i = 0; i < len; i++) 
	    Rprintf("%u ", probe[i]);
	Rprintf("\n");
    }
    
    for (i = 0; i < len; i++) {
	inner_match = 0;
	
	for (j = 0; j < current_node->current_children; j++) {
	    /** proceed to the next number in the key. **/
	    if (current_node->children[j].name == probe[i]) {
		current_node = (current_node->children + j);
		inner_match = 1; 
		break;
	    }
	}
	
	/** if there was no match then we return null to the caller. **/
	if (!inner_match) {
	    break;
	} 
	
    }

    if (inner_match == 1 && i == len) {
	return (current_node);
    }
    else {
	return NULL;
    }
}

void add_node (unsigned long* probe, int len, double rss, struct p_node* tree) {
    if (get_node(probe, len, tree) != NULL) {
	Rprintf("error in tree, trying to add a node which already exists\n");
    }
    
    int i = 0, j = 0,k=0,inner_match = 0, outer_match = 0;
    p_node_t* current_node = tree;
    
    if (__DEBUG_TREE__) {
	Rprintf("in add_node (len = %d): ", len);
	for (i = 0; i < len; i++) 
	    Rprintf("%u ", probe[i]);
	Rprintf("\n");
    }

    for (i = 0; i < len; i++) {
	inner_match = 0;
	
	for (j = 0; j < current_node->current_children; j++) {
	    /** proceed to the next number in the key. **/
	    if (current_node->children[j].name == probe[i]) {
		current_node = (current_node->children + j);
		inner_match = 1; 
		break;
	    }
	}

	/* at this point we wish to add the node. */
	if (inner_match == 0) {
	    for (j = i; j < len; j++) {
		/** allocate more space. **/
		if (current_node->current_children >= current_node->max_children) {
		    /** memory management issue. **/
		    p_node_t* new_children = Calloc(2*current_node->max_children, p_node_t);
		    for (k = 0; k < current_node->current_children; k++) {
			new_children[k] = current_node->children[k];
		    }
		    current_node->max_children = 2*current_node->max_children;
		    Free(current_node->children);
		    current_node->children = new_children;
		}
		
		/** add. **/
		/*		p_node_t* new_node = make_new_node(probe[j], ((j == (len - 1)) ? rss : -1)); */

		p_node_t* new_node = (current_node->children + current_node->current_children++);
		new_node->children = (p_node_t*) Calloc(BASE_SIZE, p_node_t);
		new_node->current_children = 0;
		new_node->max_children = BASE_SIZE;
		new_node->name = probe[j];
		new_node->rss = ((j == (len - 1)) ? rss : -1);


		current_node = new_node;

	    }
	    break;
	}
    }
}

void delete_tree(p_node_t* tree) {
    if (tree == NULL) {
	/** noop */
    }
    else if (tree->children == NULL || tree->current_children == 0) {
	Free(tree->children);
	tree->children = NULL;
    }
    else {
	int i = 0; 
	for (i = 0; i < tree->current_children; i++) 
	    delete_tree(&tree->children[i]);
	Free(tree->children);
        tree->children = NULL;
    }
}


