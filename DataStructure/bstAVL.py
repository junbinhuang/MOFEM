# BST -- AVL Tree
import sys, os
sys.path.append(os.path.dirname(sys.path[0])) 

class new_node:
    def __init__(self,key):
        self.key=key
        self.parent=None
        self.left=None
        self.right=None
        self.height=0
        self.skew=0
    
    def __str__(self):
        if self.parent and self.left and self.right:
            return '[Key, Parent, Left, Right] = '+'['+str(self.key)+', '+\
                str(self.parent.key)+', '+str(self.left.key)+', '+str(self.right.key)+']'+'\n'+\
                    '[Height, Skew] = '+'['+str(self.height)+', '+str(self.skew)+']'
        elif self.parent and self.left:
            return '[Key, Parent, Left, Right] = '+'['+str(self.key)+', '+\
                str(self.parent.key)+', '+str(self.left.key)+', '+'None'+']'+'\n'+\
                    '[Height, Skew] = '+'['+str(self.height)+', '+str(self.skew)+']'
        elif self.parent and self.right:
            return '[Key, Parent, Left, Right] = '+'['+str(self.key)+', '+\
                str(self.parent.key)+', '+'None'+', '+str(self.right.key)+']'+'\n'+\
                    '[Height, Skew] = '+'['+str(self.height)+', '+str(self.skew)+']'
        elif self.left and self.right:
            return '[Key, Parent, Left, Right] = '+'['+str(self.key)+', '+\
                'None'+', '+str(self.left.key)+', '+str(self.right.key)+']'+'\n'+\
                    '[Height, Skew] = '+'['+str(self.height)+', '+str(self.skew)+']'
        elif self.parent:
            return '[Key, Parent, Left, Right] = '+'['+str(self.key)+', '+\
                str(self.parent.key)+', '+'None'+', '+'None'+']'+'\n'+\
                    '[Height, Skew] = '+'['+str(self.height)+', '+str(self.skew)+']'
        elif self.left:
            return '[Key, Parent, Left, Right] = '+'['+str(self.key)+', '+\
                'None'+', '+str(self.left.key)+', '+'None'+']'+'\n'+\
                    '[Height, Skew] = '+'['+str(self.height)+', '+str(self.skew)+']'
        elif self.right:
            return '[Key, Parent, Left, Right] = '+'['+str(self.key)+', '+\
                'None'+', '+'None'+', '+str(self.right.key)+']'+'\n'+\
                    '[Height, Skew] = '+'['+str(self.height)+', '+str(self.skew)+']'
        else:
            return '[Key, Parent, Left, Right] = '+'['+str(self.key)+', '+\
                'None'+', '+'None'+', '+'None'+']'+'\n'+\
                    '[Height, Skew] = '+'['+str(self.height)+', '+str(self.skew)+']'

def find(node,k):
    if node.key==k:
        return node
    elif k<node.key and node.left:
        return find(node.left,k)
    elif k>node.key and node.right:
        return find(node.right,k)
    return None

def find_min(node):
    while node.left: node=node.left
    return node

def find_max(node):
    while node.right: node=node.right
    return node

def traversal(node): # In-order
    """Generator function that returns the sorted list of all nodes."""
    if node.left: yield from traversal(node.left)
    yield node
    if node.right: yield from traversal(node.right)

def is_right_child(node):
    return node.parent and (node.parent.right is node)

def is_left_child(node):
    return node.parent and (node.parent.left is node)

def successor(node):
    if node.right: return find_min(node.right)
    while is_right_child(node): node=node.parent
    return node.parent

def predecessor(node):
    if node.left: return find_max(node.left)
    while is_left_child(node): node=node.parent
    return node.parent

def update(node):
    left_height=node.left.height if node.left else -1
    right_height=node.right.height if node.right else -1

    node.height=1+max(left_height,right_height)
    node.skew=right_height-left_height

def maintain(node):
    update(node)
    balance(node)
    if node.parent:
        maintain(node.parent)

def balance(node):
    if node.skew==2:
        if node.right.skew==-1:
            right_rotate(node.right)
        left_rotate(node)
    if node.skew==-2:
        if node.left.skew==1:
            left_rotate(node.left)
        right_rotate(node)

def left_rotate(node):
    tempNode=new_node(node.key)

    tempNode.left=node.left
    if node.left: node.left.parent=tempNode
    tempNode.right=node.right.left
    if node.right.left: node.right.left.parent=tempNode
    tempNode.parent=node

    node.key=node.right.key
    if node.right.right: node.right.right.parent=node
    node.right=node.right.right
    node.left=tempNode

    update(tempNode)
    update(node)

def right_rotate(node):
    tempNode=new_node(node.key)

    tempNode.left=node.left.right
    if node.left.right: node.left.right.parent=tempNode
    tempNode.right=node.right
    if node.right: node.right.parent=tempNode
    tempNode.parent=node

    node.key=node.left.key
    node.right=tempNode
    if node.left.left: node.left.left.parent=node
    node.left=node.left.left

    update(tempNode)
    update(node)

def insert(node,k): # Empty tree is allowed.
    if node.key: # The key can be just a number so "is not None" is explicitly given.
        if k<=node.key:
            if node.left:
                insert(node.left,k)
            else:
                node.left=new_node(k)
                node.left.parent=node
                maintain(node)
        else:
            if node.right:
                insert(node.right,k)
            else:
                node.right=new_node(k)
                node.right.parent=node
                maintain(node)
    else:
        node.key=k

def delete(node):
    if node.left and node.right:
        succ=find_min(node.right)
        node.key=succ.key
        delete(succ)
    elif node.left: 
        node.key=node.left.key
        node.right=node.left.right
        node.left=node.left.left
        
        if node.left: node.left.parent=node
        if node.right: node.right.parent=node
        maintain(node)
    elif node.right:
        node.key=node.right.key
        node.left=node.right.left
        node.right=node.right.right

        if node.left: node.left.parent=node
        if node.right: node.right.parent=node
        maintain(node)
    else:
        if node.parent and is_left_child(node):
            node.parent.left=None
            maintain(node.parent)
        elif node.parent and is_right_child(node) :
            node.parent.right=None
            maintain(node.parent)
        else:
            node.key=None # Delete the root.
            maintain(node)

if __name__=='__main__':
    pass