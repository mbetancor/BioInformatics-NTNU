class SuffixTree:
    # This is a Node class that is internal to the suffixTree class.

    class Node:
        def __init__(self , value):
            self.value = value
            self.children = []
            self.isSuffixNode = True

    def __init__(self):
        self.root = self.Node("")
        self.suffixNodes = [self.root]

    def updateTree(self, string):
        for character in string:
            newSuffixNodes = []
            for suffixNode in self.suffixNodes:
                hit = next((child for child in suffixNode.children if child.value == character), None)
                if hit:
                    newSuffixNodes.append(hit)
                    hit.isSuffixNode = True
                else:
                    newNode = self.Node(character)
                    suffixNode.children.append(newNode)
                    newSuffixNodes.append(newNode)
                suffixNode.isSuffixNode = False
                    
            self.suffixNodes[1:] = newSuffixNodes

    def isInTree(self, string):
        suffixesFound = []
        output = ""
        pointer = self.root
        for character in string:
            hit = next((child for child in pointer.children if child.value == character), None)
            if hit:
                output+=character
                pointer = hit
                if hit.isSuffixNode:
                    suffixesFound.append(output)
            else:
                break

        return suffixesFound

        
# def main():
#     s = input("Enter a list of numbers: ")
#     lst = s.split()
    
#     tree = suffixTree()
    
#     for x in lst:
#         tree.insert(float(x))
        
#     for x in tree:
#         print(x)

# if __name__ == "__main__":
#     main()