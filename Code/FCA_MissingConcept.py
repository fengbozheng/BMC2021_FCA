#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 14:19:59 2021

@author: fengbozheng
"""
###############################################################################
#Step1: Mappings
###############################################################################
import csv
replaceMap = {} #key will be replaced by value while performing FCA on sequences
replaceMapReverse = {} #key could be changed to values while validating sequences
#e.g., while peforming FCA, tumor will be mapped to neoplasm
#e.g., while validating, neoplasm will be replaced by 
#[tumor, other synonym, tumor before noromalization and other synonym before normalization]
file1 = open("SynonymReplacementMapping.txt","r")
for lines1 in file1:
    line1 = lines1.split("\n")[0].split("\t")
    if replaceMap.get(line1[0],"Default") == "Default":
        replaceMap[line1[0]] = line1[1]
    else:
        print("error")
    if replaceMapReverse.get(line1[1],"Default") == "Default":
        replaceMapReverse[line1[1]] = [line1[0]]
    else:
        replaceMapReverse[line1[1]].append(line1[0])
file1.close()
#
file2 = open("NormalizationReplacementMapping.csv","r")
reader2 = csv.reader(file2)
for row2 in reader2:
    if replaceMapReverse.get(row2[0],"DEFAULT") == "DEFAULT":
        replaceMapReverse[row2[0]] = row2[1:]
    else:
        for items in row2[1:]:
            if items not in replaceMapReverse.get(row2[0]):
                replaceMapReverse[row2[0]].append(items)
file2.close()
############################################################################### 
#Step2 Intersection & Associative Role Comparison
############################################################################### 
import csv
import string
allowed = string.ascii_letters + string.digits
def check_naive(myString):
    return any(c in allowed for c in myString)
conceptLexical = {}
file1 = open("conceptInformation_1908_normalized.txt","r",encoding="ISO-8859-1")
for lines1 in file1:
    line1 = lines1.split("\n")[0].split("\t")
    lexicalItems = line1[2].lower().split("()()")
    qualified = []
    for item1 in lexicalItems:
        if check_naive(item1) != False:
            qualified.append(item1)
    #qualifiedIter = qualified[:]
    for i in range(len(qualified)):
        if replaceMap.get(qualified[i],"DEFAULT") != "DEFAULT":
            qualified[i] = replaceMap.get(qualified[i])
    conceptLexical[line1[0]] = qualified
file1.close()
#
import networkx as nx  
file22 = open("hierarchicalRelation(ParentChild)_1908.txt","rb")
Dirgraph2 = nx.read_edgelist(file22, create_using = nx.DiGraph(),nodetype = str)
def findSubhierarchy(node):
    Dic = dict(nx.bfs_successors(Dirgraph2,node))
    d = Dic.values()
    c = []
    for e in d:
        c.extend(e)
    return c 
file22.close()
###############################################################################
def lcs(S,T):
    m = len(S)
    n = len(T)
    counter = [[0]*(n+1) for x in range(m+1)]
    longest = 0
    lcs_set = set()
    for i in range(m):
        for j in range(n):
            if S[i] == T[j]:
                c = counter[i][j] + 1
                counter[i+1][j+1] = c
                if c > longest:
                    lcs_set = set()
                    longest = c
                    lcs_set.add(tuple(S[i-c+1:i+1]))
                elif c == longest:
                    lcs_set.add(tuple(S[i-c+1:i+1]))
    return lcs_set
###############################################################################
def seq_intersect(S,T):
    resultSet = set()
    interSeq = []
    for i in range(len(S)):
        if S[i] in T:
            interSeq.append(S[i])
    interSeq2 = []
    for j in range(len(T)):
        if T[j] in S:
            interSeq2.append(T[j])
    if len(interSeq)>0:
        resultSet.add(tuple(interSeq))
    if len(interSeq2)>0:
        resultSet.add(tuple(interSeq2))
    return resultSet
#
def subList(a,b): #check if a in b
    j = 0
    for i in range(len(b)):
        if b[i:i+len(a)] == a:
            j = j+1
    if j ==0:
        return False
    else:
        return True    
############################################################################### 
import csv
#import itertools
import string
allowed = string.ascii_letters + string.digits
def check_naive(myString):
    return any(c in allowed for c in myString)
#
allNameInNCIt = {}
NCItLabel = {}
file2 = open("conceptInformation_1908_normalized.txt","r",encoding="ISO-8859-1")
for lines2 in file2:
    line2 = lines2.split("\n")[0].split("\t")
    BagOfWord = line2[2].lower().split("()()")
    checkedBOW1 = []
    for item1 in BagOfWord:
        if check_naive(item1) != False:
            checkedBOW1.append(item1)
    for i in range(len(checkedBOW1)):
        if replaceMap.get(checkedBOW1[i],"DEFAULT") != "DEFAULT":
            checkedBOW1[i] = replaceMap.get(checkedBOW1[i])
    checkedBOW = tuple(checkedBOW1)
    if allNameInNCIt.get(checkedBOW,"Default") == "Default":
        allNameInNCIt[checkedBOW] = [line2[0]]
    else:
        allNameInNCIt[checkedBOW].append(line2[0])
    NCItLabel[line2[0]] = line2[3]
file2.close()
#
import networkx as nx
file2 = open("hierarchicalRelation(ChildParent)_1908.txt","rb")
Dirgraph1 = nx.read_edgelist(file2, create_using = nx.DiGraph(), nodetype = str)
file2.close()
#
def findAncestors(node):
    Dic = dict(nx.bfs_successors(Dirgraph1,node))
    d = Dic.values()
    c = []
    for e in d:
        c.extend(e)
    return c 
allNodesInGraph = set(Dirgraph1.nodes)
#
associatedRoles = {}
file1 = open("conceptRelationInferred_1908State.txt","r")
for lines1 in file1:
    line1 = lines1.split("\n")[0].split("\t")
    if line1[1]!= "ISA2019FZ":
        if associatedRoles.get(line1[0],"DEFAULT") == "DEFAULT":
            associatedRoles[line1[0]] = [(line1[1],line1[2])]
        elif (line1[1],line1[2]) not in associatedRoles.get(line1[0]):
            associatedRoles[line1[0]].append((line1[1],line1[2]))     
file1.close()
#
def checkSingleGeneral(attributeB, attributeA):  #check if B is more general than A
    if attributeB[0] in allNodesInGraph and attributeA[0] in allNodesInGraph:
        if (attributeB[0] == attributeA[0]) or (attributeB[0] in findAncestors(attributeA[0])):
            if (attributeB[1] in findAncestors(attributeA[1])) or (attributeB[1] == attributeA[1]):
                return True
            else:
                return False
        else:
            return False
    else:
        if (attributeB[0] == attributeA[0]):
            if (attributeB[1] in findAncestors(attributeA[1])) or (attributeB[1] == attributeA[1]):
                return True
            else:
                return False
        else:
            return False
#       
def roleIntersection(attributeLB, attributeLA):
    intersection = []
    for singleAttributeB in attributeLB:
        if any(checkSingleGeneral(singleAttributeB,x) for x in attributeLA):
            intersection.append(singleAttributeB)
    for singleAttributeA in attributeLA:
        if any(checkSingleGeneral(singleAttributeA,y) for y in attributeLB):
            intersection.append(singleAttributeB)
    return(list(set(intersection)))
############################################################################### 
#Step3: Detecting Missing Concepts
############################################################################### 
roots = list(Dirgraph2.successors("C2991")) 
for root in roots:
    print(root)
    allConcept = findSubhierarchy(root)
    allConcept.append(root)
    originalSet = set() #all sequences 
    #allBOW = set() #all bag of words
    for eachConcept in allConcept:
        itsSeq = conceptLexical.get(eachConcept)
        #itsBOW = tuple(sorted(set(itsSeq)))
        originalSet.add(tuple(itsSeq))
        #allBOW.add(itsBOW)       
    #
    initialSet = originalSet.copy()  # used to store all formalzied concepts
    oldIteration= set()
    newlyAdded = originalSet.copy()  # used to store newly added concepts for each iteration
    flag = 0
    while flag ==0:
        #currentIteration = initialSet.copy()
        lastIteration = initialSet.copy()     #concept names before an iteration
        oldIterationL = list(oldIteration)    #keep old concept names (no need to do intersection for them)
        newlyAddedL = list(newlyAdded)
        #last iteration new vs last iteration new
        for i in range(0,len(newlyAddedL)-1):
            for j in range(i+1,len(newlyAddedL)):
                a = lcs(newlyAddedL[i],newlyAddedL[j])
                for newlyFormed in a:
                    initialSet.add(newlyFormed)
        #last iteration old vs last iteration new        
        for k in range(0,len(newlyAddedL)):
            for l in range(0,len(oldIterationL)):
                b = lcs(newlyAddedL[k],oldIterationL[l])
                for newlyFormed in b:
                    initialSet.add(newlyFormed)
        #update variables
        newlyAdded = initialSet - lastIteration
        oldIteration = lastIteration.copy()     #keep old concept names
        if len(newlyAdded) == 0:
            flag = 1
    #
    allNewlyAdded = initialSet - originalSet
    #
    filename = NCItLabel.get(root)+"_NewConcepts_LCS(>0).csv"
    output5 = open(filename,"w")
    writer5 = csv.writer(output5)
    for eachNewOne in allNewlyAdded:
        directSubtypes = []
        directSuperTypes = []
        for existConcept in allConcept:
            existConceptLex = conceptLexical.get(existConcept)
            #if set(existConceptLex).issuperset(set(eachNewOne)):
            if subList(list(eachNewOne),existConceptLex):
                directSubtypes.append(existConcept)
            if subList(existConceptLex,list(eachNewOne)):
                directSuperTypes.append(existConcept)
        allSubtypes = "\t".join([x+": "+NCItLabel.get(x) for x in directSubtypes])
        lowerBound = directSubtypes[:]
        for child1 in directSubtypes:
            if any(x in findAncestors(child1) for x in directSubtypes):
                lowerBound.remove(child1)
        lowerBoundL = "\t".join([y+": "+NCItLabel.get(y) for y in lowerBound])
        lowerBoundRole = [associatedRoles.get(x,[]) for x in lowerBound]
        intersectionFlag = 0
        initialIntersect = lowerBoundRole[0]
        for i in range(1,len(lowerBoundRole)):
            initialIntersect = roleIntersection(initialIntersect,lowerBoundRole[i])
        if len(initialIntersect)>0:
            intersectionFlag = 1
        #
        upperBound = directSuperTypes[:]
        for parent1 in directSuperTypes:
            if any(parent1 in findAncestors(x) for x in directSuperTypes):
                upperBound.remove(parent1)
        upperBoundL = "\t".join([y+": "+NCItLabel.get(y) for y in upperBound])        
        writer5.writerow((lowerBoundL,upperBoundL,len(lowerBound),len(upperBound),intersectionFlag,allSubtypes)+eachNewOne)
    output5.close()   
