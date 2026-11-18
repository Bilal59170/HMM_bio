# Try to clear the context using GPT

### Observations $(y_j)_j$: present-day DNA variants

- **DNA is a sequence of bases**—A, C, G and T. At a given position (*locus*), different chromosome copies may carry different bases. These alternatives are called **alleles**.
- Each person has **two copies of each autosomal chromosome**, one inherited from each parent. Thus, four diploid individuals provide **$S=8$ haploid sequences** (chromosome copies).
- We compare these copies at the **same genomic position**, not at different dates. We observe present-day alleles, **not the ancient mutation events directly**.

### Notation in the paper

The observations are represented as

$$
y=\{(s_j,A_j)\}_{j=1}^{k},
$$

Symbol & Meaning:
- **$|y|=k$** is the **number of observed variant sites**, not a site itself. Their genomic positions are $s_1,\ldots,s_k$.
- **$S$** is the **number of sampled DNA sequences (chromosome copies)** studied.
- **$A_j\in\{0,1\}^{S}$** gives the alleles across those $S$ copies at variant site $s_j$.

Some remarks:
- $s_j \in [1, L)$ with $L$ the length of the DNA, but $(s_j) \subsetneq [1, L).$ The variant site are observed and taken, we don't have a variant at each position. 
- One can picture the allele data as an **$S\times k$ matrix**: rows are chromosome copies, columns are variant sites.
- It's binary because some sites are **biallelic site**, i.e. there are only two alternatives (A and G e.g). -0 and 1 label the two alleles. For example,

$$
(A,G,A,A)\quad\longrightarrow\quad A_j=(0,1,0,0).
$$

These labels are site-specific: **1 does not always mean G**.

Thus, the intuitive observation **$Y_s$** is the allele pattern across sampled copies at genomic position $s$; the paper records these patterns as **$A_j$** at variant positions **$s_j$**. Observed stretches without variants also provide information.

### Understanding the latent variable, here the trees

**Nodes / Leafs / Branches**
- **Leafs:** The leafs are the $\mathcal{S}$ observed DNA sequences.  
- **Nodes:** An internal node represents a common ancestor of the sampled DNA copies—more precisely, a coalescence event, where two ancestral lineages merge when traced backward. **Its numerical value is its age**, usually measured in generations before the present. 
- **Branch:** a branch’s length represents the elapsed evolutionary time between its child and parent nodes: $$branch_{length}=t_{parent}−t_{child}.$$
For example: 
```
             Ancestor, age 100
               /          \
              /            Ancestor, age 40
             /                /      \
          A, age 0         B, age 0  C, age 0
```

**Link with mutations: some example**
Imagine tracing each copy backward through its ancestors. When two ancestral paths reach the same ancestral DNA copy, they join. This is a **coalescence event**, represented by an internal node.

Here is an illustrative tree: 
```
                      R  (1,000 generations ago)
                     / \
                    /   \
     P (200 generations) Q (400 generations)
             / \               / \
            /   \             /   \
           1     2           3     4
                    Present: u = 0
```
- Copies 1 and 2 share an ancestral copy 200 generations ago.
- Copies 3 and 4 share one 400 generations ago.
- All four share an ancestral copy 1,000 generations ago.

Suppose at a variant site we observe
$$A_j=(1,1,0,0).$$
One possible explanation is a mutation on the branch between $R$ and $P$: copies 1 and 2 inherited that change, while copies 3 and 4 did not. Mutations occur along branches; internal nodes represent shared ancestry, not mutations. 

However, several tree can explain a given observation. Here, for example this other tree:

```
                     R (1,000 generations ago)
                    / \
                   /   \
 P (200 generations)     \
          / \             \
         1   2             ├── 3
                           └── 4

Present: u = 0
```

-> Hence it's needed to observe several locuses to get more accurate trees. 

**Trees and position $s$ in DNA sequence**
- We observe several DNA sequences, and each position $s$ has a hidden local tree describing their ancestry.
- A tree exists even at positions with no observed variation. If $A_s$ records the DNA variant in each sampled copy, a variable position might look like $(0,1,0,1)$; an invariant position might be $(0,0,0,0)$.
- We do not have one tree per mutation as the same tree can apply across a stretch of DNA. The 2 key concepts are 
	- **Mutations**: occur on branches and provide evidence about the underlying tree. 
	- **Recombination**: can change the underlying tree as we move along the genome. 
Hence, a local tree explain the mutations while the change of tree along the locuses explains the recombination. 

Here is some scheme of the tree. 
```
Genomic position → ──────────────────────────────

Local genealogy:   [ Tree A       ][ Tree B       ]
Observed mutations:    *    *          *   *   *
```
-> Let's remark that a given tree, e.g tree A, spans locuses with mutations ( $*$ symbol) and without mutation (no $*$)


### The smoothing problem

For a hidden state $X_s,$ the $n_p$​ particles $X_s^{(1)},…,X_s^{(n_p)}$​ are candidate realizations used to approximate its distribution. Here, we are interested in a ==smoothing problem== because we want to infer trees across the genome, not only the tree at the current position: 
$$p(X_{0:s} | Y_{0:s}).$$
Smoothing estimates $p(X_u | Y_{0:s})$ for $u < s$ at (biological) time $s.$ The particles of the earlier tree $X_y$ are updated using observations collected **after** that position. For example, a mutation observed at position 200 can help us infer the tree at position 100, because the two positions may share the same tree—or closely related genealogies.

**Some comparison with the robot case:**

| Robot tracking                                                                                             | Genomic inference                                                             |
| ---------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------- |
| Index $t$: clock time                                                                                      | Index $s$: position along the DNA                                             |
| Hidden state $X_t$​: robot position                                                                        | Hidden state $X_s$: local genealogical tree                                   |
| Observation $Y_t$: sensor reading                                                                          | Observation $Y_s$: genetic data at that position                              |
| Motion model                                                                                               | Model of how recombination changes the genealogy                              |
| Smoothing at time $t$: recover the whole trajectory of the robot $X_{0:t}$ observing the sensors $Y_{0:t}$ | Smoothing at time $s$: recover all the tree up until the genome of index $s$. |
