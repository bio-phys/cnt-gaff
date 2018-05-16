with open('out.pdb','r') as a:
    with open('system.pdb','w') as w:
        prev = '1'
        for line in a:
            l = line.split()
            if l[5] != prev:
                w.write("TER\n")
            w.write(line)
            prev = l[5]
