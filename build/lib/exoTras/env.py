Adata = None
Inter_adata = None
Thershold = None
Same = None

def initializer(inter_adataI, sameI, thersholdI):
    global Inter_adata
    global Thershold
    global Same
    Inter_adata = inter_adataI
    Thershold = thersholdI
    Same = sameI

def initializer_simple(inter_adataI):
    global Inter_adata
    Inter_adata = inter_adataI

def initializer_adata(adataI, sameI, thersholdI):
    global Adata
    global Thershold
    global Same
    Adata = adataI
    Thershold = thersholdI
    Same = sameI

