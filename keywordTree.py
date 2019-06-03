class suffixTree:
	# This is a Node class that is internal to the suffixTree class.
	
	def __init__(self):
		self.suffixNodes = [0,1]

	class Node:
		def __init__(self,data):
			self.data = data
			self.children = []
			self.isSuffixNode = True



		def addChildren(self, obj):
			self.children.append(obj)




class KeyWordTree:
	# This is a Node class that is internal to the BinarySearchTree class. 
	class Node:
		def __init__(self,val):
			self.val = val
			self.out = {}

			
		def getVal(self):
			return self.val
		
		def setVal(self,newval):
			self.val = newval
			
		def getLeft(self):
			return self.left
		
		def getRight(self):
			return self.right
		
		def setLeft(self,newleft):
			self.left = newleft
			
		def setRight(self,newright):
			self.right = newright
			
		# This method deserves a little explanation. It does an inorder traversal
		# of the nodes of the tree yielding all the values. In this way, we get
		# the values in ascending order.
		def __iter__(self):
			if self.left != None:
				for elem in self.left:
					yield elem
					
			yield self.val
			
			if self.right != None:
				for elem in self.right:
					yield elem

		def __repr__(self):
			return "BinarySearchTree.Node(" + repr(self.val) + "," + repr(self.left) + "," + repr(self.right) + ")"            
			
	# Below are the methods of the BinarySearchTree class. 
	def __init__(self, root=None):
		self.root = root
		
	def insert(self,val):
		self.root = BinarySearchTree.__insert(self.root,val)
		
	def __insert(root,val):
		if root == None:
			return BinarySearchTree.Node(val)
		
		if val < root.getVal():
			root.setLeft(BinarySearchTree.__insert(root.getLeft(),val))
		else:
			root.setRight(BinarySearchTree.__insert(root.getRight(),val))
			
		return root
		
	def __iter__(self):
		if self.root != None:
			return iter(self.root)
		else:
			return iter([])

	def __str__(self):
		return "BinarySearchTree(" + repr(self.root) + ")"
 
def main():
	s = input("Enter a list of numbers: ")
	lst = s.split()
	
	tree = BinarySearchTree()
	
	for x in lst:
		tree.insert(float(x))
		
	for x in tree:
		print(x)

if __name__ == "__main__":
	main()