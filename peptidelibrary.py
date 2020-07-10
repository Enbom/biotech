from operator import add

#[C, H, N, O, S]
aminoAcids = {"A":[3, 7, 1, 2, 0], "R":[6, 14, 4, 2, 0], "N":[4, 8, 2, 3, 0], "D":[4, 7, 1, 4, 0], "C":[3, 7, 1, 2, 1], "Q":[5, 10, 2, 3, 0], "E":[5, 9, 1, 4, 0], "G":[2, 5, 1, 2, 0], "H":[6, 9, 3, 2, 0], "I":[6, 13, 1, 2, 0], "L":[6, 13, 1, 2, 0], "K":[6, 14, 2, 2, 0], "M":[5, 11, 1, 2, 1], "F":[9, 11, 1, 2, 0], "P":[5, 9, 1, 2, 0], "S":[3, 7, 1, 3, 0], "T":[4, 9, 1, 3, 0], "W":[11, 12, 2, 2, 0], "Y":[9, 11, 1, 3, 0], "V":[5, 11, 1, 2, 0]   }

def elementalMassPercentage(seq, element):
    masses=peptideElementalAverageMass(seq)
    elementalComposition = {"C":0,"H":0,"N":0,"O":0,"S":0}
    for i in range(0,5):
        elementalComposition["CHNOS"[i]]=masses[i]/sum(masses)
    return elementalComposition

    if element == "C":
        return masses[0]/sum(masses)
    elif element == "H":
        return masses[1]/sum(masses)
    elif element == "N":
        return masses[2]/sum(masses)
    elif element == "O":
        return masses[3]/sum(masses)
    elif element == "S":
        return masses[4]/sum(masses)
    else:
        return False

def peptideElementalAverageMass(seq): #atomic units >>
    composition=peptideComposition(seq)
    return [12.011* composition[0],1.008* composition[1],14.007* composition[2],15.99* composition[3],32.06* composition[4]]

def peptideComposition(seq):
    composition = [0, -(len(seq)-1)*2, 0, -(len(seq)-1), 0] #Adjusting for peptide bonds, -H2O per peptide bond
    for char in seq:
        composition = map(add, composition, aminoAcidElementalComposition(char) )
    return list(composition) #[C, H, N, O, S]


def aminoAcidElementalComposition(aminoAcid):
    return aminoAcids[aminoAcid] #[C, H, N, O, S]


print(elementalMassPercentage("ELPKLPDDKVLIRSRSNCPKGKVWNGFDCKSPFAFS", "C"))
