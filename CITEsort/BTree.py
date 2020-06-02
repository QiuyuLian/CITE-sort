#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 18 15:41:14 2019

@author: lqyair
"""

class BTree:

    def __init__(self, key, left = None, right = None, indices = None, stop=None,\
                 all_clustering_dic = None, where_dominant = None,weight=None,ll=None,bic=None):
        self.key = key # a str
        self.right = right # a BstNode
        self.left = left # a BstNode
        self.indices = indices # a []
        self.all_clustering_dic = all_clustering_dic
        self.weight = weight
        self.ll = ll
        self.bic = bic
        #self.marker_summary = marker_summary # a pd.df
        #self.para = para # a {} parameters for qualified components
        self.where_dominant = where_dominant # str ("left"/"right"), indicator of edge color
        self.stop = stop # legacy
            
    
    def display(self):
        lines, _, _, _ = self._display_aux()
        for line in lines:
            print(line)

    def _display_aux(self):
        """Returns list of strings, width, height, and horizontal coordinate of the root."""
        # No child.
        if self.right is None and self.left is None:
        #if self.right.key is 'leaf' and self.left.key is 'leaf':
            line = '%s' % '_'.join(self.key)
            width = len(line)
            height = 1
            middle = width // 2
            return [line], width, height, middle

        # Only left child.
        if self.right is None:
        #if self.right.key is 'leaf':
            lines, n, p, x = self.left._display_aux()
            s = '%s' % '_'.join(self.key)
            u = len(s)
            first_line = (x + 1) * ' ' + (n - x - 1) * '_' + s
            second_line = x * ' ' + '/' + (n - x - 1 + u) * ' '
            shifted_lines = [line + u * ' ' for line in lines]
            return [first_line, second_line] + shifted_lines, n + u, p + 2, n + u // 2

        # Only right child.
        if self.left is None:
        #if self.left.key is 'leaf':
            lines, n, p, x = self.right._display_aux()
            s = '%s' % '_'.join(self.key)
            u = len(s)
            first_line = s + x * '_' + (n - x) * ' '
            second_line = (u + x) * ' ' + '\\' + (n - x - 1) * ' '
            shifted_lines = [u * ' ' + line for line in lines]
            return [first_line, second_line] + shifted_lines, n + u, p + 2, u // 2

        # Two children.
        left, n, p, x = self.left._display_aux()
        right, m, q, y = self.right._display_aux()
        s = '%s' % '_'.join(self.key)
        u = len(s)
        first_line = (x + 1) * ' ' + (n - x - 1) * '_' + s + y * '_' + (m - y) * ' '
        second_line = x * ' ' + '/' + (n - x - 1 + u + y) * ' ' + '\\' + (m - y - 1) * ' '
        if p < q:
            left += [n * ' '] * (q - p)
        elif q < p:
            right += [m * ' '] * (p - q)
        zipped_lines = zip(left, right)
        lines = [first_line, second_line] + [a + u * ' ' + b for a, b in zipped_lines]
        return lines, n + m + u, max(p, q) + 2, n + u // 2


    
    

