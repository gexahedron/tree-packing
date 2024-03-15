CXX ?= g++
CXXFLAGS += -std=c++11 -O3
INCLUDES += -I .
# PROGRAM_NAME = graceful_tree_packing
# PROGRAM_NAME = mod_graceful_tree_packing
# PROGRAM_NAME = planar_tree_packing
# PROGRAM_NAME = all_planar_tree_packing
CC = g++
FLAGS = -Wall -O3 -std=c++11 -DNDEBUG -I.
FLOBJ = $(FLAGS) -c
FLEXE = $(FLAGS) -o $@ $^
TIE_BREAK_POLICY = neighbor_future_degree_tie_break_policy.o tie_break_policy.o
TRIANGULATION = $(TIE_BREAK_POLICY) triangulation.o tie_break_min_fill_triangulation.o
TREE_DECOMPOSITION = tree_decomposition.o tie_break_min_fill_triangulation_tree_decomposition.o triangulation_tree_decomposition.o
STUFF = child_node.o cluster_node.o cluster_tree.o light_graph.o order.o timer.o vertex_set.o $(TRIANGULATION) $(TREE_DECOMPOSITION)
# PROGRAM_NAME = all_series_parallel_graphs_packing
PROGRAM_NAME = all_bipartite_graphs_packing

all :
	${CXX} $(CXXFLAGS) $(INCLUDES) $(PROGRAM_NAME).cpp $(FLEXE) -o $(PROGRAM_NAME)

#all_series_parallel_graphs_packing: child_node.o cluster_node.o cluster_tree.o light_graph.o order.o  timer.o vertex_set.o $(TRIANGULATION) $(TREE_DECOMPOSITION) all_series_parallel_graphs_packing.cpp
#	$(CC) $(FLEXE)

# all_sp_graphs_no_t2_no_t1_packing: child_node.o cluster_node.o cluster_tree.o light_graph.o order.o  timer.o vertex_set.o $(TRIANGULATION) $(TREE_DECOMPOSITION) all_sp_graphs_no_t2_no_t1_packing.cpp
# 	$(CC) $(FLEXE)

# all_sp_graphs_no_stars_packing: child_node.o cluster_node.o cluster_tree.o light_graph.o order.o  timer.o vertex_set.o $(TRIANGULATION) $(TREE_DECOMPOSITION) all_sp_graphs_no_stars_packing.cpp
# 	$(CC) $(FLEXE)

#all_leaf_graphs_packing: child_node.o cluster_node.o cluster_tree.o light_graph.o order.o  timer.o vertex_set.o $(TRIANGULATION) $(TREE_DECOMPOSITION) all_leaf_graphs_packing.cpp
#	$(CC) $(FLEXE)

#%.o: %.cpp %.h
#	$(CC) $(FLOBJ) $<

