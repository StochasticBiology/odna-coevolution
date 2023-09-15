#### next q -- can we add a new node to the most advanced possible existing one, so the branches don't all come from zero?
# next q -- transitions embedded. idea -- first attach source nodes then dest nodes on subgraph rooted at appropriate source
# q -- this looks like it's working, but what if we have the same destination observation that needs to be reached from two different sources? -- "awkward" -- currently throws error
# q -- TB drug data creates loops / multiple parents? -- seems to be fixed by taking unique observations.
# q -- working on general case. 1,1,0,0,1,0,1,0,0,1 throws an issue -- cautiously optimistic?
 
import csv
import sys

# is dest compatible with src, in the sense that we can in principle have a src -> dest path?
def compatible(src, dest):
    for i in range(len(src)):
        if src[i] == '1' and dest[i] == '0':
            return False
    return True

# initialise lists for destination and source observations
data = []
sourcedata = []

# based on command line, do we have just a destination file or a source file too?
if(len(sys.argv) > 2):
    sourcefile = 1
else:
    sourcefile = 0

# at least open and read the destination file
fp = open(sys.argv[1], "r")

for line in fp.readlines():
    tmpline = line.rstrip("\n").replace(",","")
    data.append(tmpline)

# either read the source file, or give us an empty list of sources
if sourcefile == 1:
    fp = open(sys.argv[2], "r")
    for line in fp.readlines():
        tmpline = line.rstrip("\n").replace(",","")
        sourcedata.append(tmpline)
else:
    sourcedata = []

# identify any observations that are only present in the source, not the destination, set
sourceonlys = [x for x in sourcedata if x not in data]

# root the tree at 0^L, and add to dataset
root = "0"*len(data[1])
data.append(root)
sourcedata.append(root)

# get the required source/ancestor states for each observation. if we don't have a source file, set all of these to root (independent, cross-sectional observations)
sources = {}
if sourcefile == 1:
    for i in range(0,len(data)):
        if data[i] in sources:
            print("I've already got a source ("+sources[data[i]][0]+") for "+data[i]+" -- this will break the code at the moment")
            sources[data[i]].append(sourcedata[i])
        else:
            sources[data[i]] = [sourcedata[i]]
else:
    data = list(set(data))
    for i in range(0,len(data)):
        sources[data[i]] = [root]

# set required source for any source-only observations to be the root itself
for item in sourceonlys:
    sources[item] = [root]

# construct our to-do list of observations to attach
todo = list(data+sourceonlys)

# initialise branch count and edge set
branch = 0
edges = []

print("Sources: "+str(sources))

# while we still have observations to include
while len([x for x in todo if x != root]) > 0:
    print("\nBranch "+str(branch))
    print("TO DO: "+str(todo))
    branch = branch+1
    prevadd = "none"
    # initialise compatibility dictionary -- indexed list of compatible sets of remaining observations
    compat = {}
    compatall = {}
    # loop through pairs of remaining observations assessing compatibility
    for dp in todo:
        compat[dp] = []
        for dp2 in todo:
            if dp != dp2 and compatible(dp2, dp):
                compat[dp].append(dp2)
    # initialise loop
    loopstate = 0
    # get indices for compatibility index -- ie remaining observations
    newset = [x for x in compat]
    # loop through the next longest available "necklace" 
    while loopstate != 2:
        print("TO DO: "+str(todo))
        # get number of compatible observations for each remaining observation
        ncomp = [len(compat[x]) for x in newset]
        # find the observation with the most compatible observations
        try:
            maxi = ncomp.index(max(ncomp))
        except:
            print("Error here. newset is "+str(newset)+"\nncomp is "+str(ncomp))
            exit()
        # flag this for adding, and consider only the observations compatible with this one in the next iteration
        newadd = [x for x in newset][maxi]
        print(" -- going to add "+newadd)
        # is this a destination observation? if so, grab its source node
        # to allow a dest node to have multiple sources, we process the required sources one at a time, first popping the first remaining one
        required = sources[newadd][0]
        if newadd != root:
            sources[newadd].pop(0)
        newset = compat[newadd]
        newokset = []
        print("newset is "+str(newset)+" bc compat[newadd] is "+str(compat[newadd]))
        # impose compatibility with source node
        for member in newset:
            if compatible(required, member):
                print(member+" is compatible with required "+required)
                newokset.append(member)
            else:
                print(member+" is not compatible with required "+required)
        newset = newokset
        print("  Choosing "+newadd+" with "+str(max(ncomp))+" compatible strings ("+str(len(newset))+" after requirement of ", str(required)+")")
        # if newset is empty, then either we have a required source that doesn't appear in the data, or our source is already in the DAG
        # (without(?)) any additional states that can be threaded between the source and our current observation. so we need to just
        # attach to that source
        addedexisting = 0
        if len(newset) == 0 and newadd != root:
            possibles = [x[1] for x in edges]
            if required in possibles:
                print("Found "+required+" in existing DAG, adding "+str([required,newadd])+" and "+str([newadd,prevadd]))
                # we risk adding required and newadd but not connecting to prevadd; BUT if prevadd is root (how can this happen?) we get a different bug
                # I think this happens if this module is met on the first trip through a branch. trying to catch with an empty prevadd marker
                if prevadd != "none":
                    edges.append([newadd,prevadd])
                edges.append([required,newadd])
                addedexisting = 1
                loopstate = 2
            else:
                print("Error: Can't find "+required)
                exit()
                
        # if the next to be added would be the root, see if we attach later to an existing path
        if newadd == root:
            bestscore = 0
            bestnew = root
            possibles = [x[1] for x in edges]
            print("    possibles:"+str(possibles))
            for possnew in possibles:
                if compatible(possnew, prevadd) and possnew.count("1") > bestscore:
                    bestscore = possnew.count("1")
                    bestnew = possnew
            print("    chose "+bestnew)
            newadd = bestnew
            loopstate = 2
        else:
            todo.remove(newadd)
        # add to growing list            
        # if this is a leaf, update loop state; otherwise, add an edge to the next lowest observation
        if loopstate == 0:
            loopstate = 1
        elif addedexisting == 0:
            edges.append([newadd, prevadd])
            print("Adding "+str([newadd, prevadd]))
        prevadd = newadd

# write edge list to CSV file
with open(sys.argv[1]+"-outs2.csv", "w", newline="") as f:
    f.write("From,To\n")
    writer = csv.writer(f)
    writer.writerows(edges)
