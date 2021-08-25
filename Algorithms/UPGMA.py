import numpy as np

def get_shortest_dist(inp):
	current_min = inp[0, 1]
	current_min_coords = (0, 1)
	for i in range(inp.shape[0] - 1):
		for j in range(i + 1, inp.shape[1]):
			if inp[i, j] < current_min:
				current_min = inp[i, j]
				current_min_coords = (i, j)
	return current_min, current_min_coords

def calculate_dist(dist_1, dist_2):
	res = (dist_1[1] * dist_1[0]) + (dist_2[1] * dist_2[0])
	return res / (dist_1[1] + dist_2[1])

def merge_nodes(nodes, dist):
	new_el = (nodes[0]['id'], nodes[1]['id'])
	new_node = create_node(new_el, dist, nodes[0], nodes[1])
	new_node['size'] = new_node['left']['size'] + new_node['right']['size']
	if new_node['left']['size'] > 1:
		new_node['dist_left'] = (dist / 2) - (new_node['left']['dist'] / 2)
	if new_node['right']['size'] > 1:
		new_node['dist_right'] = (dist / 2) - (new_node['right']['dist'] / 2)
	return new_node 

def merge(inp, clusters, els, dist):
	n_clusters = [merge_nodes((clusters[els[0]], clusters[els[1]]), dist)]
	for i in range(len(clusters)):
		if i not in els:
			n_clusters.append(clusters[i])
	n_matrix = np.zeros((len(n_clusters), len(n_clusters)))
	original = [0]
	for i in range(len(inp)):
		if i not in els:
			original.append(i)
	for i in range(n_matrix.shape[0] - 1):
		for j in range(i + 1, n_matrix.shape[1]):
				if i == 0:
					n_val = calculate_dist((inp[original[j], els[0]], clusters[els[0]]['size']), (inp[original[j], els[1]], clusters[els[1]]['size']))
				else:
					n_val = inp[original[i], original[j]]
				n_matrix[i, j] = n_val
				n_matrix[j, i] = n_val
	return n_matrix, n_clusters
					

# Create a node that contains elements, distance and connections, if any
def create_node(el, dist, left=None, right=None):
	node = {}
	node['id'] = el
	node['dist'] = dist
	node['size'] = 1
	if not left is None and not right is None:
		node['left'] = left
		node['right'] = right
	return node

# Converts number to letter
# Long numbers will be represented by multiple letters
# Examples:
# 	0 -> a
#	2 -> c
#	27 -> aa
def number_to_str(num):
	res = chr(ord('a') + int(num) % 26)
	while num // 26 != 0:
		num = (num // 26) - 1
		res = chr(ord('a') + int(num) % 26) + res
	return res

def UPGMA(inp):
	# Create clusters
	clusters = [create_node(number_to_str(i), 0) for i in range(len(inp))]
	current = inp
	# While there is more than one cluster
	while len(current) >= 2:
		c_min, c_ind = get_shortest_dist(current)
		current, clusters = merge(current, clusters, c_ind, c_min)
	return clusters[0]
