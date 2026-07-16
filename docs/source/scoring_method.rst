How PRS Scoring Works
=====================

This page explains, in plain terms, how ``aoutools`` turns a weights file and
the *All of Us* genomic data into a polygenic risk score (PRS). More to the
point, it explains the handful of subtle things the library does to make sure
that number is correct. You do not need any of this to use the library; the
:doc:`how-to guides <guides/index>` are enough. Read on if you want to trust the
scores, reproduce them, or understand why a value looks the way it does.

Every behaviour described here is checked against real *All of Us* data before
each release.

Two symbols recur below. :math:`w` is a variant's **weight**, and :math:`n` is
the number of **non-reference** copies a person carries at that variant
(:math:`0`, :math:`1`, or :math:`2`).


The basic idea
--------------

A PRS is a weighted sum. A weights file lists, for each variant, an **effect
allele** and a weight (an effect size from a GWAS). Each person's score is

.. math::

   \text{PRS} \;=\; \sum_{\text{variants}} w \times c,
   \qquad
   c = \text{copies of the effect allele } (0,\ 1,\ \text{or } 2).

A person carries two copies of every position, one from each parent, so the copy
count :math:`c` is 0, 1, or 2. If :math:`w = 0.4` and someone carries the effect
allele on both chromosomes, that variant adds :math:`0.4 \times 2 = 0.8` to their
score.

The arithmetic really is that simple. The rest of this page is about one thing:
getting the copy count :math:`c` right for every person at every variant. That is
where genomic data hides several traps, each able to produce a wrong but
perfectly plausible number.


How the data is stored
----------------------

Several of the traps below turn on how *All of Us* physically stores the
genotypes, so it is worth setting out first. At each variant position, the people
who match the reference are not written out one by one. They are recorded
together as a single **reference block**, which says, in effect, "everyone not
listed here matches the reference." Only the people who carry an alternate allele
get an individual entry, and that entry lists just the alleles they actually
have. A person whose genotype could not be determined also gets an entry, marked
as a **no-call**.

A person's copy count at a variant (the :math:`c` above) is read from their
entry. Reference-matching people have no entry, so their count is not read but
inferred: they match the reference, which is 2 reference copies. A no-call person
does have an entry, but it says "unknown", so they are scored as nothing. Keeping
those two cases apart is what the whole design turns on, and Traps 2 and 4 below
are where it matters most.


Trap 1: The effect allele can be either the reference or the alternate
----------------------------------------------------------------------

Genomic data describes each variant relative to a **reference genome**. There is
a reference allele (REF, the base the reference carries) and one or more
alternate alleles (ALT, the bases that differ). A weights file is not built that
way. The GWAS that wrote it chose whichever allele it wanted as the effect
allele, with no promise that this is the ALT. For roughly half of all variants,
the effect allele is actually the reference.

``aoutools`` works this out **per variant**, by looking at the actual reference
at that position. It never assumes a whole file is oriented one way. An earlier
design offered a single file-wide switch ("the reference is the effect allele"),
but it was removed: it silently dropped every variant that disagreed with the
switch, and a normal PGS Catalog file mixes both orientations.

Which allele is the effect allele decides how the copy count is read:

.. math::

   c \;=\;
   \begin{cases}
     n, & \text{effect allele is the ALT (count non-reference copies)} \\[4pt]
     2 - n, & \text{effect allele is the REF (see Trap 2)}
   \end{cases}

**Example.** One weights file, effect weight :math:`w = 0.3`, used at two
variants. At the first, the effect allele ``A`` is the *alternate*; at the
second, the effect allele ``C`` is the *reference*.

.. math::

   \begin{array}{l l l}
   \textbf{data (REF/ALT)} & \textbf{effect allele} & \textbf{a carrier's contribution} \\
   \hline
   \text{G / A} & \text{A = ALT} & \text{G/A} \rightarrow c = 1 \rightarrow 0.3 \times 1 = 0.3 \\
   \text{C / T} & \text{C = REF} & \text{C/T} \rightarrow c = 2-1 = 1 \rightarrow 0.3 \times 1 = 0.3 \\
                &                & \text{C/C} \rightarrow c = 2-0 = 2 \rightarrow 0.3 \times 2 = 0.6 \\
   \end{array}

Both rows came from the same file. The library decided each one on its own, from
the data, with no file-wide setting involved.


Trap 2: When the effect allele is the reference, most people are invisible
--------------------------------------------------------------------------

This is the most important point — and the least obvious.

Recall from `How the data is stored`_ that a reference-matching person has no
entry at a variant. They are not stored as "0 copies"; they are simply absent, so
no per-person calculation ever visits them.

When the effect allele is the ALT, this is fine: an absent person carries zero
copies of the ALT, and zero is exactly what they should contribute. When the
effect allele is the **reference**, it is a problem. An absent person carries
:math:`2 - 0 = 2` copies of the effect allele, the maximum possible — yet there
is no record of them to add up.

``aoutools`` solves this with a small piece of algebra. For a reference-effect
variant the true contribution is :math:`w\,(2 - n)`. The library splits it into a
part that is the same for everyone and a part that depends on the person:

.. math::

   \underbrace{\,w\,(2 - n)\,}_{\text{true contribution}}
   \;=\;
   \underbrace{\,2w\,}_{\substack{\text{row-level constant} \\ \text{added to everyone}}}
   \;-\;
   \underbrace{\,w\,n\,}_{\substack{\text{per person} \\ \text{zero when there is no entry}}}

The per-person term :math:`w\,n` is already zero for an absent (reference-
matching) person, so it needs no special handling. The :math:`2w` term is the
same constant for everyone, so the library adds it once, at the variant level, to
every person. That row-level constant is the only thing that can reach people who
have no entry to visit.

**Example** with :math:`w = 0.5`:

.. math::

   \begin{array}{l c c l}
   \textbf{genotype} & n & \textbf{copies } (2-n) & \textbf{contribution } (2w - w n) \\
   \hline
   \text{ref / ref (no entry)} & 0 & 2 & 2(0.5) - 0.5(0) = 1.0 \\
   \text{heterozygous}         & 1 & 1 & 2(0.5) - 0.5(1) = 0.5 \\
   \text{homozygous alternate} & 2 & 0 & 2(0.5) - 0.5(2) = 0.0 \\
   \end{array}

The ref/ref person is stored nowhere — yet still receives the correct
:math:`1.0`.
That value comes from the shared constant added to everyone, not from a record
read for them, because there is none.

Three quiet rules keep this honest, and all three are tested:

- the constant is added only for variants actually **found** in the data. Adding
  it for an absent variant would invent signal for the whole cohort.
- it is **cancelled** for a person whose genotype at that variant is *unknown*
  (see Trap 4). An unknown genotype is not the same as a reference match.
- it is genuinely the same for everyone, so it can never reorder people.


Trap 3: Sites with more than one alternate allele
-------------------------------------------------

At some positions the population carries more than one alternate. Say the
reference is ``C`` and both ``G`` and ``T`` occur: this is a **multi-allelic**
site. About one in five *All of Us* variants is multi-allelic, so it is not a
corner case.

``aoutools`` first **splits** such a site into separate two-allele rows, one per
alternate. A ``C/G,T`` site becomes a ``C/G`` row and a ``C/T`` row, so each
alternate can be matched and scored on its own. (This is Hail's ``split_multi``
step.)

Splitting introduces a subtlety. On the row for one alternate, say ``T``, the
split rewrites the *other* alternate ``G`` as if it were the reference. A person
who actually carries ``C/G`` therefore looks, on the ``T`` row, exactly like
someone who matches the reference. That is correct when the effect allele is
``T`` (they carry no ``T``), but wrong when the effect allele is the reference
``C``, because they really carry one ``C``, not two.

To get both cases right, the library keeps **two** counts for each person at each
variant:

- the copies of *this row's specific alternate*, used when the effect allele is
  the ALT; and
- the person's **total** non-reference count :math:`n`, taken from *before* the
  split, used when the effect allele is the REF.

**Example.** Reference ``C``, alternates ``G`` and ``T``. The weight's effect
allele is the reference ``C``, with :math:`w = 0.2`. Take a person carrying
``C/G``, who has exactly one ``C`` and should contribute :math:`0.2 \times 1`:

.. math::

   \begin{array}{l l}
   \textbf{count used on the T row} & \textbf{result} \\
   \hline
   \text{this row's count (G looks like REF, so genotype looks like C/C)}
       & 0.2 \times 2 = 0.4 \quad \text{(wrong)} \\
   \text{pre-split total } (n = 1,\ \text{so } 2 - n = 1)
       & 0.2 \times 1 = 0.2 \quad \text{(correct)}
   \end{array}

The library uses the pre-split total for reference-effect weights, so the ``C/G``
carrier scores :math:`0.2`, not :math:`0.4`. The mistake is invisible at ordinary
two-allele sites, where the two counts are equal, which is exactly why it went
unnoticed until it was tested at real multi-allelic sites.


Trap 4: An unknown genotype is not a reference match
----------------------------------------------------

Sometimes the data records that a person's genotype at a variant could not be
determined. This is a **no-call**, and it is not the same as a reference match. A
reference match is *known*: the person has the common allele. A no-call is
*unknown*.

``aoutools`` treats a no-call as contributing nothing: not two reference copies,
not zero-with-an-offset. In particular, the :math:`2w` constant from Trap 2 is
explicitly removed for no-call people, so an unknown genotype is never quietly
scored as two copies of the reference. Confusing these two cases was a real bug
in an earlier version, and the code now guards against it specifically.

**Example**, reference-effect weight :math:`w = 0.5` (the same setup as Trap 2):

.. math::

   \begin{array}{l l}
   \text{ref / ref} & \text{known to carry 2 copies} \rightarrow 0.5 \times 2 = 1.0 \\
   \text{no-call}   & \text{genotype unknown} \rightarrow 0
       \quad (\text{the } 2w \text{ constant is added, then cancelled back out})
   \end{array}

Without the cancellation, the no-call person would wrongly collect the :math:`2w`
constant and look identical to someone confirmed to carry two reference copies.


Trap 5: The same variant can be written more than one way
---------------------------------------------------------

Before a weights row can be scored, it has to be matched to the right variant in
the data, and the same variant can be **spelled** differently on each side.
``ATG/A`` and ``TG`` both describe a deletion of ``TG``. A GWAS might list a SNP
as ``A/G`` while the data spells it ``AAAG/GAAG``. If matching used the raw text,
these would fail to line up and the variant would silently score nothing.

``aoutools`` avoids that in two steps.

- **Reduce to a canonical spelling.** Every variant, on both sides, is reduced to
  its **minimal representation**: shared padding letters are trimmed until the
  alleles are as short as possible. ``AAAG/GAAG`` reduces to ``A/G``, and
  ``ATG/A`` reduces to a clean ``TG`` deletion, so the two sides agree. This
  matters even after the split in Trap 3, because the split step reduces only the
  rows it actually splits. An already-two-allele variant written non-minimally,
  like ``AAAG/GAAG``, would otherwise slip through un-reduced, so ``aoutools``
  reduces every row itself.
- **Match on the unordered allele set.** Two variants match when they sit at the
  same position and carry the same pair of alleles, order ignored. This is what
  lets ``G/A`` in the file match ``A/G`` in the data (Trap 1). Only after they
  match does the library decide which allele is the effect allele. Matching on
  *ordered* alleles instead used to drop every variant whose alleles sorted the
  "wrong" way, about half of them.

A weights row that names the *same* allele as both its effect and non-effect
allele, a degenerate ``A/A`` row, describes no real variant. Its two alleles are
identical, so the set can never equal a genuine two-allele site. Such rows match
nothing and contribute zero, instead of matching something by accident.

One direction of trimming is safe; the other needs a guard. Trimming shared
letters off the **end** of both alleles keeps the variant at the same position,
and it is how a GWAS names most indels, so the library relies on it. Trimming
shared letters off the **front** would move the variant to a different position,
and genomic tools cannot do that safely: they can only raise an error or silently
discard the variant. ``aoutools`` chooses to raise a loud error rather than drop
it silently. No variant of this shape exists in the current *All of Us* data
(zero out of six million checked), so this is a tripwire for a hypothetical
future data release — not something you will hit today.


A gallery of variant shapes
---------------------------

The table below pulls the whole page together: a small catalog of variant
shapes, one per row. Each shows a variant as *All of Us* stores it, the effect
allele a weight names for it (and whether that allele is the reference or the
alternate, resolved against the stored variant as in Trap 1), and the scoring
situation the pair creates. Between them they cover every path the scoring code
takes.

.. list-table::
   :header-rows: 1
   :widths: 12 20 16 52

   * - Locus
     - Alleles (REF/ALT)
     - Effect allele
     - What it exercises
   * - ``chr1:1000``
     - ``A/G``
     - ``G`` (ALT)
     - The simplest case: a two-allele SNP whose reference base sorts before its
       alternate. Four people carry ``G``; everyone else is stored as a
       reference block.
   * - ``chr1:2000``
     - ``A/G``
     - ``A`` (REF)
     - A reference-effect SNP where one person has a **no-call**. An unknown
       genotype must score nothing, not two reference copies (Trap 4).
   * - ``chr1:3000``
     - ``A/T``
     - ``A`` (absent)
     - The data carries ``A/T`` here, but a weight names ``A/G``. That weighted
       variant is **not in the data at all**, so it must contribute nothing to
       anyone, not even a reference offset.
   * - ``chr1:4000``
     - ``A/T``
     - ``G`` (absent)
     - The same missing variant as the row above, but the weight now names the
       *other* allele as the effect allele. An absent variant must still score
       nothing whichever way round it is written.
   * - ``chr1:5000``
     - ``C/G/T``
     - ``C`` (REF)
     - A **multi-allelic** site whose people carry different alternates. After
       the site is split, a carrier of the alternate the weight did not name
       must not be miscounted (Trap 3).
   * - ``chr1:6000``
     - ``AGGGC/A/GGGGC``
     - ``G`` (ALT)
     - A multi-allelic site written non-minimally. Trimming the shared ending
       turns ``AGGGC/GGGGC`` into a plain ``A/G`` SNP at the same position, so
       the ``G``-effect weight still matches (Trap 5).
   * - ``chr1:7000``
     - ``A/C/G``
     - ``A`` (REF)
     - A multi-allelic site where one person is homozygous for ``C``, the
       alternate the weight did *not* name. This is the worst case for the
       reference offset: mishandling it credits two copies the person does not
       carry (Trap 3).
   * - ``chr1:8000``
     - ``G/A``
     - ``A`` (ALT)
     - A SNP whose reference base sorts *after* its alternate. Order-blind
       matching is what keeps it from being silently dropped (Traps 1 and 5).
   * - ``chr1:8500``
     - ``AAAG/GAAG``
     - ``G`` (ALT)
     - A two-allele variant written non-minimally. Because it already has only
       two alleles, the split step leaves it untouched, so the library has to
       reduce it to ``A/G`` itself before the ``G``-effect weight matches
       (Trap 5).
   * - ``chr1:9000``
     - ``A/C``
     - ``C`` (ALT)
     - A plain SNP carrying a negative, non-integer weight, so the total is real
       arithmetic rather than a copy count.


A note on the weights themselves
--------------------------------

The weight :math:`w` is just a number, and it can be negative: a protective
allele legitimately lowers the score. ``aoutools`` does the real arithmetic
rather than counting copies, so negative and non-integer weights are handled as
they are.

One conversion is worth knowing about. Some GWAS publish effect sizes as **odds
ratios** instead of the log-odds that a sum assumes. If you pass
``log_transform_weight=True``, the library takes the natural logarithm of each
weight first (:math:`w \rightarrow \ln w`), so the scores add up on the correct
scale. Use it only when your weights really are odds ratios.


Running efficiently on a very large dataset
-------------------------------------------

The *All of Us* genomic dataset is enormous, and reading all of it for every
score would be slow and expensive. ``aoutools`` keeps the work affordable in
three ways, none of which changes the numbers:

- **Chunking.** The weights are processed in **chunks**. For each chunk the
  library reads back only the exact positions those variants sit at, then adds
  the per-chunk results together. This is the main reason the library is
  practical on the full dataset, and the final score is identical whether you use
  one chunk or many.
- **Batch scoring.** Computing several scores with ``calculate_prs_batch`` reads
  the VDS once for all of them, rather than once per score. Each score comes out
  identical to computing it on its own, so batching is a speed optimisation, not
  a different calculation.
- **Scoring a subset of people.** If you only need a defined cohort, pass
  ``samples_to_keep`` and the library restricts the calculation to those people
  before it counts anything. Each person's PRS depends only on their own
  genotypes, so narrowing the set never changes anyone else's score.

Both the chunk equality and the batch-equals-single equality are tested.


What ``aoutools`` does **not** do
---------------------------------

Two things are deliberately left to you, because guessing at them would do more
harm than good:

- **Strand and genome build.** The library assumes your weights are on the same
  DNA strand and the same genome build (GRCh38) as the *All of Us* data. It does
  not detect or correct a strand flip.
- **Palindromic variants.** Variants whose two alleles are strand-symmetric
  (``A/T`` and ``C/G``) are not flagged. If your weights and the data disagree on
  strand for such a variant, the alleles alone cannot reveal it, so you must
  resolve it upstream.

Make sure your weights file matches the reference build and strand before
scoring.


In short
--------

You do not have to manage any of this yourself. Give ``aoutools`` a weights file
on the GRCh38 build and the correct strand, and it handles allele orientation,
multi-allelic sites, reference-matching and no-call genotypes, and alternate
spellings of the same variant, so the copy count behind every score is the right
one. The one thing left to you is making sure your weights are on the same build
and strand as the data before you start.

Each behaviour above is validated against the real *All of Us* data before every
release, so the scores you get are the scores the method intends.
