# SMC Project: Make things more clear

- [SMC Project: Make things more clear](#smc-project-make-things-more-clear)
  - [Understanding the sequential coalescent with recombination model](#understanding-the-sequential-coalescent-with-recombination-model)
    - [Vocabulary \& Position of the problem](#vocabulary--position-of-the-problem)
    - [Observations (y\_j)\_j: present-day DNA variants](#observations-y_j_j-present-day-dna-variants)
    - [Latent variable x, sequence of tree (and more)](#latent-variable-x-sequence-of-tree-and-more)
    - [Recombination](#recombination)
      - [Recombination and the exit rate](#recombination-and-the-exit-rate)
    - [Coalescence](#coalescence)
    - [Markov jump process from Poisson and Cox processes](#markov-jump-process-from-poisson-and-cox-processes)
        - [Main idea](#main-idea)
        - [Poisson process, uniform rate lambda](#poisson-process-uniform-rate-lambda)
        - [Poisson process, non-uniform rate lambda(s)](#poisson-process-non-uniform-rate-lambdas)
        - [Cox process with marks: equations (8) and (9)](#cox-process-with-marks-equations-8-and-9)
    - [Structure of pi\_x(x)](#structure-of-pi_xx)
      - [What x and x\_s are (Appendix)](#what-x-and-x_s-are-appendix)
      - [Joint model](#joint-model)
  - [The smoothing problem](#the-smoothing-problem)
  - [Filters: course \& qualitative explanation](#filters-course--qualitative-explanation)
    - [Typical algorithm, FK, and Filters](#typical-algorithm-fk-and-filters)
    - [Expression of G\_t, choice of M\_t](#expression-of-g_t-choice-of-m_t)
    - [Bootstrap vs Guided (qualitative)](#bootstrap-vs-guided-qualitative)
  - ["Methods" section of the paper](#methods-section-of-the-paper)
    - [Particle filter (algo 1)](#particle-filter-algo-1)
    - [CMJP (algo 2)](#cmjp-algo-2)
    - [Choice of waypoints (empty)](#choice-of-waypoints-empty)
    - [Lokkahead likelihood (empty)](#lokkahead-likelihood-empty)
    - [EM and Variationnal Bayes for MLE (empty)](#em-and-variationnal-bayes-for-mle-empty)
  - [Other parts](#other-parts)


## Understanding the sequential coalescent with recombination model

### Vocabulary & Position of the problem
We would like to construct genealogical trees froms DNA sequences. 

- **DNA is a sequence of bases**—A, C, G and T. At a given position (*locus*), different chromosome copies may carry different bases. These alternatives are called **alleles**.
- Each person has **two copies of each autosomal chromosome**, one inherited from each parent. Thus, four diploid individuals provide **$S=8$ haploid sequences** (chromosome copies).
- We compare these copies at the **same genomic position**, not at different dates. We observe present-day alleles, **not the ancient mutation events directly**.


### Observations (y_j)_j: present-day DNA variants


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
### Latent variable x, sequence of tree (and more)

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




### Recombination

For one autosomal chromosome pair, the father has copies \(f_1,f_2\) and the mother has \(m_1,m_2\). Each passes **one** copy to the child, which may contain segments from both of their own copies:

```text
Father:  f₁  A A A A A A      Mother:  m₁  C C C C C C
         f₂  B B B B B B               m₂  D D D D D D
Child:   c_f A A A B B B               c_m C C D D D D
```

Thus $$c_f=\operatorname{mix}(f_1,f_2) \quad c_m=\operatorname{mix}(m_1,m_2).$$ **Recombination** is the exchange of segments that can produce these patchwork copies; it is not a fusion of the parents’ chromosomes.


#### Recombination and the exit rate

Recombination change the ancester while going through the DNA. For example, going through $c_f,$ we switch from $A$ to $B.$ Therefore, **a recombination can change the ancestry we trace** as we move along the genome, so the **local tree can change**. Take the example of 3 DNA sequences. Instead of grouping copies 1 and 2 together in the first tree, we group copies 1 and 3 together in the second. 

An **exit rate** tells us how quickly a Markov process leaves its current state. Intuitively, the more recombination, the quicker a tree changes for another. Here, the state is the local tree \(x_s\), and the process moves along genomic position \(s\), not evolutionary time. Thus, \(\rho(s)B(x_s)\) is the rate of a recombination event on the current tree per nucleotide of genomic distance.

The units are:

$$
\underbrace{\rho(s)}_{\text{per nucleotide per generation}}
\times
\underbrace{B(x_s)}_{\text{total branch length in generations}} =
\underbrace{\rho(s)B(x_s)}_{\text{per nucleotide}}.
$$

Over a very short genomic distance \(\Delta s\), the probability of a recombination event is approximately

$$
\Pr\!\left(
\text{recombination event in } [s,s+\Delta s)
\mid x_s
\right)
\approx
\rho(s)B(x_s)\,\Delta s.
$$

A tree with more total branch length offers more ancestral opportunities for recombination. 

### Coalescence

Look at the genealogy **backwards in time**. Any two people (or, more precisely, two copies of a gene) have a most recent common ancestor. Going back from the present, their two ancestral lineages eventually **merge** into one. That merger is a **coalescence**, and the time at which it happens is the **coalescence time**.

```
present:   A     B     C
           |     |     |
           |_____|     |     <- A and B coalesce (recent)
              |        |
              |________|     <- AB and C coalesce (older)
                  |
              common ancestor of all three
```

A genealogical tree is a sequence of such mergers. In the paper, "evolutionary time" $u$ is the number of generations before the present.

**What is $N_e$?**

$N_e$ is the **effective population size**. It is the size of an idealised population (random mating, constant size, no selection) that would show the same amount of genetic drift as the real population. It is usually smaller than the census size. For example, humans have an $N_e$ of roughly $10{,}000$ even though billions of people are alive today. The number reflects the population's long-term history, including bottlenecks.

For a diploid species, a population of size $N_e$ has $2N_e$ gene copies at each locus.

**Why is the coalescence rate $C(u) = 1/(2N_e(u))$?**

Take two gene copies in the same generation. Each one descends from a copy in the previous generation. The chance that they pick the **same** parent copy is $1/(2N_e)$.That chance is the probability that the two lineages coalesce in one generation, so it is the **pairwise coalescence rate per generation**.
- With $b$ lineages (ancêtre/lignée en français) present, there are many possible pairs. The rate at which a specific lineage merges with any of the other $b-1$ is therefore about $(b-1)\,C$.
- A small population has a high $C$, because lineages find each other quickly.
- A large population has a low $C$, because lineages stay apart for a long time.

For constant $N_e$, two lineages take on average $2N_e$ generations to coalesce. With $N_e = 10{,}000$ that is about $20{,}000$ generations, or $C = 5\times 10^{-5}$.

**Why does $C$ depend on $u$?**

Population size changes over time (expansions, bottlenecks), so $N_e(u)$ is a function of the generation $u$ in the past. This is why the paper's inference of demography amounts to estimating the function $C(u)$. Once you have it, you recover $
N_e(u) = \frac{1}{2C(u)}.$

### Markov jump process from Poisson and Cox processes

##### Main idea
A Markov jump process stays in a state, here a tree, for a random length of the genome. Then it jumps to a new state. The density of a whole path is built from two pieces:

- Staying: the probability of no jump during a stretch, which is $\exp{-\text{rate}\times\text{length}}$.
- Jumping: the density that a jump happens at a given place, which is the rate at that place.

For example, for a constant rate $q$ the waiting time until the next jump has density $q,e^{-qt}$. The formula for $\pi_x(x)$ is this idea repeated for every jump.
##### Poisson process, uniform rate lambda
A Poisson process with rate $\lambda>0$ on an interval $[0,L)$ is a random set of points (events) such that:

1. Independent increments: the numbers of events in disjoint intervals are independent.
2. Poisson counts: the number of events in an interval of length $\ell$ follows a Poisson law,
$$
P(N=k) = e^{-\lambda \ell}\,\frac{(\lambda \ell)^k}{k!}
$$
3. No simultaneous events: events occur one at a time.


Using that given $N=k,$ the positions are i.i.d. uniform, with density $1/\ell$ each, we can deduce the probability of ocurrence of sequence $(s_1 , \dots , s_k):$

$$
p(s_1,\dots,s_k) = \underbrace{e^{-\lambda \ell}\frac{(\lambda \ell)^k}{k!}}_{\text{count proba}} \cdot \underbrace{\ell^{-k}}_{\text{positions proba}} = \frac{1}{k!}\,e^{-\lambda \ell}\,\lambda^k.
$$

##### Poisson process, non-uniform rate lambda(s)
Here, $\lambda(\cdot)$ is not constant over $[1, L).$ Hence, the mean value and count probability are respectively
$$
\Lambda = \int_1^L \lambda(s)\,ds, \qquad P(N=k) = e^{-\Lambda}\,\frac{\Lambda^k}{k!}.
$$

We just replace $\ell L$ by $\Lambda.$ The same way as before we have th probability of ocurrence. Given $N=k$ are i.i.d. with density $\lambda(s)/\Lambda,$

$$
p(s_1,\dots,s_k) = e^{-\Lambda}\frac{\Lambda^k}{k!}\dot \quad \prod_{i=1}^k \frac{\lambda(s_i)}{\Lambda} = \frac{1}{k!}\,e^{-\Lambda}\prod_{i=1}^k \lambda(s_i)
$$

The $\Lambda^k$ cancels, and $\lambda \to \lambda(s_i)$ because each point sits at a specific location.

##### Cox process with marks: equations (8) and (9)

**Setup.** Each event is a pair $(\tilde s_i, m_i)$: a locus $\tilde s_i \in [1,L)$ and a label (mark) $m_i \in M$, where $M$ is a discrete set. In the paper, $m$ is the allele pattern $A \in \{0,1\}^S \setminus \{(0,\dots,0)\}$ observed in the $S$ samples.
$$
y = \big((\tilde s_1, m_1), \dots, (\tilde s_k, m_k)\big), \qquad |y| = k
$$
**Rate per mark.** Given the hidden state $X_s = x_s$, each mark $m$ has its own rate:
$$
r(s, m \mid x_s)
$$
**Eq. (9), total rate.** This is the intensity of "an event of any type":
$$
r(x_s) \coloneqq \sum_{m \in M} r(s, m \mid x_s)
$$
For mutations:
$$
r(s, A \mid x_s) = P(A \mid x_s, \mu(s)), \qquad r(x_s) = \sum_{A \neq (0,\dots,0)} P(A \mid x_s, \mu(s)) = \mu(s)\,B(x_s)
$$
**Eq. (8), density of $y$ given the path $x$.**
$$
r(y \mid x)
= \frac{1}{|y|!}\,
\underbrace{\exp\left\{-\int_1^L r(x_s)\,ds\right\}}_{\text{no event of any type elsewhere}}\;
\underbrace{\prod_{i=1}^{|y|} r(\tilde s_i, m_i \mid x_{\tilde s_i})}_{\text{event of type } m_i \text{ at } \tilde s_i}
$$
Difference with the classic Poisson process

|             | Classic Poisson                | Marked (paper)                                                                        |
| ----------- | ------------------------------ | ------------------------------------------------------------------------------------- |
| Event       | position $s_i$ only            | position and label $(\tilde s_i, m_i)$                                                |
| Exponential | $\exp\{-\int \lambda(s)\,ds\}$ | $\exp\{-\int r(x_s)\,ds\}$, **total** rate over all $m$ (eq. 9)                       |
| Product     | $\prod_i \lambda(s_i)$         | $\prod_i r(\tilde s_i, m_i \mid x_{\tilde s_i})$, rate of the **observed** label only |

- **Exponential:** it uses the sum over all labels, because "no event" means no event of any type, including labels that are never observed.
- **Product:** it uses only the rate of the observed label $m_i$, because an event with a specific label $m_i$ happened there.
- **Reduction:** with a single label, $|M|=1$, the formulas reduce to the classic Poisson case.


### Structure of pi_x(x)

$$
\pi_x(x) =
\underbrace{\exp\left\{-\int_1^L \rho(s)\,B(x_s)\,ds\right\}}_{\text{Block 1: no other recombination}}
\prod_{j=1}^{|x|}
\Big[
\underbrace{\rho(s_j)}_{\text{2a}}\;
\underbrace{\exp\left\{-\int_{\nu_j}^{\tau_j} b_u(x_{s_j})\,C(u)\,du\right\}}_{\text{2b}}\;
\underbrace{C(\tau_j)}_{\text{2c}}
\Big]
\underbrace{(ds)^{|x|}(du)^{2|x|}}_{\text{Block 3: measure}}
$$

**Principle.** A Markov jump process is described by "survive, then jump". The density of a path is the product of the survival probabilities and the jump rates:

$$
q\,e^{-qt} = \underbrace{e^{-qt}}_{\text{survive}}\cdot\underbrace{q}_{\text{jump}}
$$

#### What x and x_s are (Appendix)

- **$x$ is more than a sequence of trees.** It also records the number of recombinations $|x|$, their loci $s_j$, and their recombination and coalescence times $\nu_j$ and $\tau_j$. These are exactly the coordinates in Block 3. The implementation keeps them to compute the sufficient statistics for inferring the coalescence and recombination rates.
- **Some recombinations do not change the tree.** $|x|$ counts them anyway. A "jump" of the process is therefore not always a visible change of local tree.
- **$x_s$ includes non-local branches.** The state is the local genealogy plus the non-local branches used by SCRM. $b_u(x_s)$ and $B(x_s)$ are computed on this full object.
- **Exact versus pruned.** The formula is exact for the CwR process if all non-local branches are kept. SCRM prunes some of them, and the formula is then applied to the pruned $x$, so it describes the approximate process actually simulated.
- **Infinite root branch.** At a recombination point $s_j$, $x_{s_j}$ is the genealogy just to the **left** of $s_j$. It includes an infinite branch above the root, so $b_u(x_{s_j}) = 1$ for $u > \text{root}(x_{s_j})$.
- **$B$ stops at the root.** $B(x_s)=\int_0^{\text{root}(x_s)} b_u(x_s)\,du$ is finite and does not include the infinite branch.

**Block 1: no recombination elsewhere (survival along the genome)**

$$
\exp\left\{-\int_1^L \rho(s)\,B(x_s)\,ds\right\}
$$

- $\rho(s)B(x_s)$ is the exit rate of the state $x_s$. A longer tree has more branch length that can be hit, so it changes more often.
- $x_s$ is piecewise constant, so the integral is a sum over segments: $\sum_{\text{states}} \rho\,B(\text{state})\times\text{segment length}$.
- $B$ is the total length of the finite branches, including the non-local branches kept in the state. The infinite root branch is excluded.

**Block 2: one factor per recombination $j=1,\dots,|x|$**

| Factor | Meaning | Type |
|---|---|---|
| $\rho(s_j)$ | Recombination at locus $s_j$, at time $\nu_j$ on the tree | jump rate |
| $\exp\{-\int_{\nu_j}^{\tau_j} b_u(x_{s_j})\,C(u)\,du\}$ | The detached lineage joins no branch between $\nu_j$ and $\tau_j$ | survival |
| $C(\tau_j)$ | It joins a given branch exactly at $\tau_j$ | jump rate |

- **2a.** My reading is that $\nu_j$ is uniform on the total branch length, with density $1/B$. This cancels the factor $B$ of the exit rate, so only $\rho(s_j)$ remains.
- **2b.** At evolutionary time $u$ there are $b_u(x_{s_j})$ branches to join, each at rate $C(u)=1/2N_e(u)$. Here $x_{s_j}$ is the genealogy **just to the left** of $s_j$.
- **Above the root.** For $u$ above the root, $b_u=1$ because of the infinite branch. The detached lineage can therefore still coalesce at $\tau_j >\text{root}(x_{s_j})$, at rate $C(u)$. In that case the lineage joins the branch above the root, and the root of the tree moves up.
- **2c.** There is no factor $b_{\tau_j}$. I read this as the choice of branch being discrete and part of the state. **The Appendix does not address this point**, so it remains an interpretation.
- **2b and 2c together** form the waiting-time density of the coalescence: survive until $\tau_j$, then coalesce.
- **No visible change.** Some of these events leave the local tree unchanged, for example a recombination on a non-local branch. They are still counted in $|x|$ and still contribute their factor.

**Block 3: reference measure**

$$
(ds)^{|x|}(du)^{2|x|}
$$

- $\pi_x$ is a **density**. Each recombination has three continuous coordinates $(s_j,\nu_j,\tau_j)$, which are encoded in $x$. This gives $|x|$ factors $ds$ and $2|x|$ factors $du$.
- The dimension depends on $|x|$, so the path space is a union of spaces of different dimensions.


#### Joint model

$$
\pi(x,y) = \pi_x(x)\; r(y\mid x), \qquad \pi(x\mid y)\propto \pi_x(x)\; r(y\mid x)
$$


## The smoothing problem

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

## Filters: course & qualitative explanation
### Typical algorithm, FK, and Filters
Typical algorithm: 
$$
\begin{array}{rl}
\hline
& \textbf{Algorithm 1:} \text{Generic PF algorithm} \\
& \text{Operations involving index } n \text{ must be performed} \\
& \text{for } n=1,\ldots,N. \\
\hline
1: & X_0^n\sim \mathbb{M}_0(dx_0) \\
2: & w_0^n\gets G_0(X_0^n) \\
3: & W_0^n\gets w_0^n\big/\sum_{m=1}^{N}w_0^m \\
4: & \textbf{for } t=1 \textbf{ to } T \textbf{ do} \\
5: & \quad A_t^{1:N}\sim\operatorname{resample}(W_{t-1}^{1:N}) \\
6: & \quad X_t^n\sim M_t(X_{t-1}^{A_t^n},dx_t) \\
7: & \quad w_t^n\gets G_t(X_{t-1}^{A_t^n},X_t^n) \\
8: & \quad W_t^n\gets w_t^n\big/\sum_{m=1}^{N}w_t^m \\
9: & \textbf{end for} \\
\hline
\end{array}
$$

The algorithm depends on four type of quantities: $\mathbb{M}_0, M_t, G_t, G_0.$ $\mathbb{M}_0, M_t$ are chosen by the user, while $G_t, G_0$ are determined by the Feynman-Kac model (cf some detail below)
FK Bootstrap model:
$$
\begin{align*}
  \mathbb{M}_0 &= \mathbb{P}_0, \quad M_t(x_{t-1}, dx_t) = P_t(x_{t-1}, dx_t), \\
  G_t(x_{t-1}, x_t) &= f_t(y_t | x_t), \quad G_0(x_0) = f_0(y_0 | x_0)
\end{align*}
$$
FK guided model:
$$
\begin{align*}
  \mathbb{M}_0(dx_0)G_0(x_0) &= \mathbb{P}_0(dx_0)f_0(y_0 | x_0), \\
  M_t(x_{t-1}, dx_t)G_t(x_{t-1}, x_t) &= P_t(x_{t-1}, dx_t)f_t(y_t | x_t).
\end{align*}
$$

### Expression of G_t, choice of M_t
You aren't supposed to work out all four. The system is underdetermined on purpose, and the free choice is exactly what the guided filter lets you design.

1. **You choose the proposals** $\mathbb{M}_0$ and $M_t$. These are design choices.
2. **The potentials then follow** from the two equations:
$$
G_0(x_0) = f_0(y_0\mid x_0)\,\frac{d\mathbb{P}_0}{d\mathbb{M}_0}(x_0),
\qquad
G_t(x_{t-1},x_t) = f_t(y_t\mid x_t)\,\frac{dP_t(x_{t-1},\cdot)}{dM_t(x_{t-1},\cdot)}(x_t).
$$
With densities, these become $G_0 = f_0\, \dfrac{p_0}{m_0}$ and $G_t = f_t\, \dfrac{p_t(x_t\mid x_{t-1})}{m_t(x_t\mid x_{t-1})}$.

The two equations only make sure that $\mathbb{M}_t \prod G_s = \mathbb{P}_t \prod f_s$. That means every valid choice targets **the same** $\mathbb{Q}_t$ and the same $L_t = p(y_{0:t})$. The choice changes the variance of the estimates, not what is being estimated.

**Requirements on the choice of $M_t$**
- **Support:** $M_t(x_{t-1},\cdot)$ must give positive probability wherever $P_t(x_{t-1},dx_t)f_t(y_t\mid x_t)$ does. Otherwise the Radon–Nikodym derivative doesn't exist.
- **Practical:** you need to be able to **sample** from $M_t$ and **evaluate** $G_t$ pointwise.

**Some possible choices**

| Choice of $M_t$ | Resulting $G_t$ |
|---|---|
| $P_t$ | $f_t(y_t\mid x_t)$: the **bootstrap** filter, a special case of guided |
| $\propto P_t(x_{t-1},dx_t)\,f_t(y_t\mid x_t)$ (optimal) | $p(y_t\mid x_{t-1})$, which doesn't depend on $x_t$ |
| An approximation of the optimal one (e.g. from a Gaussian or linearised model) | $f_t \cdot p_t/m_t$ |

So the bootstrap filter is just one solution of the guided equations, the one where you pick $M_t = P_t$.

### Bootstrap vs Guided (qualitative)
The two filters differ in **when the observation $y_t$ is used**:
- The **bootstrap** filter moves particles first and only uses $y_t$ afterwards, to score them.
- The **guided** filter uses $y_t$ while it moves them.

**Bootstrap: move blindly, then judge**

- **Move:** $X_t^n \sim P_t(X_{t-1}^{A_t^n}, dx_t)$. Particles follow the prior dynamics and ignore $y_t$.
- **Weight:** $G_t = f_t(y_t \mid X_t^n)$. Afterwards, each particle is scored by how well it explains $y_t$.

It's like shooting with your eyes closed and then checking which shots hit the target.

**When it fails:** when $f_t(y_t\mid \cdot)$ is much more **concentrated** than $P_t(x_{t-1},\cdot)$. That happens when the observation is very informative or the dynamics are very noisy. Then:
- most particles land where $f_t$ is nearly zero, so their weights are close to 0;
- a handful of particles carry almost all the weight, so the ESS collapses;
- after resampling, you end up with many copies of the same few particles, so the approximation of $\mathbb{Q}_t$ is poor.

The same failure happens when $y_t$ is an **outlier** with respect to the prior: almost no particle lands near it.

**Guided: aim using $y_t$, then correct**

- **Move:** $X_t^n \sim M_t(X_{t-1}^{A_t^n}, dx_t)$, where $M_t$ already takes $y_t$ into account. It aims at the region where $P_t f_t$ is large.
- **Weight:** $G_t = f_t \cdot \dfrac{dP_t}{dM_t}$. This corrects for the fact that particles were not drawn from $P_t$.

Because particles are already placed in the right region, the potentials $G_t$ are much **flatter** across particles. Flat weights mean a high ESS, less resampling, and more diversity among particles.

**The ideal case** is the optimal proposal $M_t(x_{t-1},dx_t) \propto P_t(x_{t-1},dx_t)\, f_t(y_t\mid x_t)$. Then $G_t(x_{t-1},x_t) = p(y_t \mid x_{t-1})$ no longer depends on the new $X_t$. The move step adds no weight variance at all.

**Why bootstrap could be interesting**

The bootstrap filter isn't simply worse. It has real advantages:
- **Simplicity:** you only need to **simulate** from $P_t$, not evaluate its density.
- **No design choice:** you don't have to build an $M_t$, and a poorly chosen $M_t$ can be *worse* than $P_t$, especially if its tails are lighter.
- **It works well when each $y_t$ is weakly informative** compared with the dynamics.

The last point applies to your paper. Each mutation carries little information about the tree, so the bootstrap-style approach is reasonable step by step. The real problem there is that the relevant information is spread **far ahead** along the genome. That's why the authors don't build a guided proposal: they add a lookahead, which corresponds to the auxiliary FK model in your notes.

## "Methods" section of the paper
### Particle filter (algo 1)

### CMJP (algo 2)
### Choice of waypoints (empty) 
### Lokkahead likelihood (empty)
### EM and Variationnal Bayes for MLE (empty)
## Other parts 

- "Particle filters": Algo of smoothing for discrete process
- "Continuous time & space Markov process": Adaptation of the algo to continuous time. 
	- Use of the RN theorem
	- Choice of the waypoints is important. Not necessarily the one observed.
- "Lookahead likelihood": For less variance, consider a kind of Auxiliary Particle Filter. 
	- Uses future information: $\tilde{p}(y_{i+1}|x_i)$
	- Uses a simplified likelihood of the observations $y$ knowing $x.$ See Apendix for details
- Choosing waypoints 
	- strategy to avoid having the $ESS$ to be too low. 
- Figure 1: ?? 
- Parameters inference :
	- $\pi_x(x)$ is seen as a continuous Markov Jump process of parameter $\theta.$
	- $\theta$ can be written as a fraction of 2 $\mathbb{E}.$ This quantity can be approximated by an EM algorithm. 
	- An improvement (in which sense?) can be made with a Variation Bayes method, where $\theta_{VB}$ is sampled from a $\Gamma$ law. 