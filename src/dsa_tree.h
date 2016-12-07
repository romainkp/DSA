
#ifndef __DSA_TREE__

#define __DSA_TREE__


/** the default size of the children vector. **/
#define BASE_SIZE 8


/** this is the prime node. each node represents a term, a path is a model **/
struct p_node {
  struct p_node* children; /*points to the children on the node*/
  int current_children; /*indicates which slot of the pointer is next to be filled*/
  int max_children; /*indicates the potential number of children that can be stored */
  unsigned long name; /* the term represented by the node */
  double rss; /* RSS for a model ending with this term */
};
typedef struct p_node p_node_t;

p_node_t* make_new_node(unsigned long name, double rss);
p_node_t* get_node (unsigned long* probe, int len, p_node_t* tree);
void add_node (unsigned long* probe, int len, double rss, struct p_node* tree);
void delete_tree(p_node_t* tree); 

#endif
