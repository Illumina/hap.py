Example where left-shifting can lead to trouble
===============================================

Consider the following REF/ALT sequence pairs.

```
REF   CTGTG------TGTGTGTGTGAAAA
ALT1  CTGTGTGTGTGTGTGTGTGTGAAAA
ALT2  CTGTGTGTGAGTGTGTGTGTGAAAA
```

=> Left-shifting ALT1

```
REF   C------TGTGTGTGTGTGTGAAAA
ALT1  CTGTGTGTGTGTGTGTGTGTGAAAA
```

... leads to this variant.

```
C -> CTGTGTG at pos. 1
```

=> Left-shifting ALT2

```
REF   CTGT------GTGTGTGTGTGAAAA
ALT2  CTGTGTGTGAGTGTGTGTGTGAAAA
```

... leads to this variant:

```
T -> TGTGTGA at pos. 4
```

This is problematic since in VCF notation, these two variants will result in two different 
records. When using `bcftools merge`, instead of one insertion we suddenly have two -- all because 
of a single-base error.

Moreover, when comparing VCF variants row by row between ALT1 and ALT2, 
we will not pick up that we actually have two mismatching variants rather than
two independent insertions (we lose the information that ALT1 and ALT2 will be never 
present together). 
