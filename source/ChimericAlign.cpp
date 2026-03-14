#include "ChimericAlign.h"

ChimericAlign::ChimericAlign(ChimericSegment &seg1in, ChimericSegment &seg2in, int chimScoreIn, Genome &genomeIn, ReadAlign *RAin)
                              : seg1(seg1in), seg2(seg2in),chimScore(chimScoreIn), RA(RAin), P(seg1in.P), pCh(P.pCh), mapGen(genomeIn) {
    stitchingDone=false;

    al1=&seg1.align;
    al2=&seg2.align;

    if (al1->roStart > al2->roStart)
        swap (al1,al2);

    // Find ex1/ex2: the exon pair that overlaps in read-space but maps to different chromosomes.
    // This identifies where the chimeric junction lives within multi-mate alignments.
    auto chrIdx = [this](uint gPos) -> uint {
        // binary search for chromosome containing gPos
        uint lo = 0, hi = mapGen.nChrReal;
        while (lo + 1 < hi) {
            uint mid = (lo + hi) / 2;
            if (mapGen.chrStart[mid] <= gPos)
                lo = mid;
            else
                hi = mid;
        }
        return lo;
    };

    int bestEx1 = -1, bestEx2 = -1;
    uint bestOverlap = 0;
    for (uint i1 = 0; i1 < al1->nExons; i1++) {
        uint rs1 = al1->exons[i1][EX_R];
        uint re1 = rs1 + al1->exons[i1][EX_L];
        uint gs1 = al1->exons[i1][EX_G];
        for (uint i2 = 0; i2 < al2->nExons; i2++) {
            uint rs2 = al2->exons[i2][EX_R];
            uint re2 = rs2 + al2->exons[i2][EX_L];
            uint gs2 = al2->exons[i2][EX_G];
            uint s = max(rs1, rs2);
            uint e = min(re1, re2);
            if (s >= e) continue; // no read-space overlap
            if (gs1 - rs1 == gs2 - rs2) continue; // same genomic mapping, not a junction
            if (chrIdx(gs1) == chrIdx(gs2)) continue; // same chromosome - splice, not chimeric
            uint overlap = e - s;
            if (overlap > bestOverlap) {
                bestOverlap = overlap;
                bestEx1 = i1;
                bestEx2 = i2;
            }
        }
    }

    // If the overlapping cross-chromosome exons are on the same mate AND both alignments
    // span multiple mates, prefer inter-mate junction (routes into "mates bracket" path).
    bool useOverlap = (bestEx1 >= 0);

    if (useOverlap && al1->exons[bestEx1][EX_iFrag] == al2->exons[bestEx2][EX_iFrag]) {
        // Same-mate overlap found. If the overlap is small (e.g. just shared scaffold),
        // and both alignments span multiple mates on different chromosomes, prefer
        // inter-mate junction. If the overlap is substantial, the within-mate junction is real.
        bool al1_multi = false, al2_multi = false;
        uint frag0 = al1->exons[0][EX_iFrag];
        for (uint i = 1; i < al1->nExons; i++)
            if (al1->exons[i][EX_iFrag] != frag0) { al1_multi = true; break; }
        frag0 = al2->exons[0][EX_iFrag];
        for (uint i = 1; i < al2->nExons; i++)
            if (al2->exons[i][EX_iFrag] != frag0) { al2_multi = true; break; }
        if (al1_multi && al2_multi && bestOverlap < pCh.junctionOverhangMin + 20)
            useOverlap = false; // small same-mate overlap + multi-mate → try inter-mate
    }

    if (useOverlap) {
        ex1 = bestEx1;
        ex2 = bestEx2;
        // Ensure iFrag ordering for chimericCheck: need iFrag1 <= iFrag2
        if (al1->exons[ex1][EX_iFrag] > al2->exons[ex2][EX_iFrag]) {
            swap(al1, al2);
            swap(ex1, ex2);
        }
    } else {
        // Look for inter-mate junction: the two alignments map to different chromosomes
        // and both span multiple mates. Determine which alignment is better on each mate
        // and place the junction at the mate boundary.
        // Check if al1 and al2 map to different chromosomes (use first exon of each).
        bool diffChrom = (al1->nExons > 0 && al2->nExons > 0 &&
                          chrIdx(al1->exons[0][EX_G]) != chrIdx(al2->exons[0][EX_G]));

        if (diffChrom) {
            // Determine which alignment is better on each iFrag by summing EX_L per mate.
            int al1_frag0 = 0, al1_frag1 = 0, al2_frag0 = 0, al2_frag1 = 0;
            for (uint i = 0; i < al1->nExons; i++) {
                if (al1->exons[i][EX_iFrag] == 0) al1_frag0 += al1->exons[i][EX_L];
                else                               al1_frag1 += al1->exons[i][EX_L];
            }
            for (uint i = 0; i < al2->nExons; i++) {
                if (al2->exons[i][EX_iFrag] == 0) al2_frag0 += al2->exons[i][EX_L];
                else                               al2_frag1 += al2->exons[i][EX_L];
            }

            // Assign: al1 should own the mate where it's relatively better, al2 the other.
            // We need iFrag(ex1) < iFrag(ex2) for the "mates bracket" path.
            // If al1 is better on frag0 and al2 on frag1: al1 keeps frag0, al2 keeps frag1 → OK
            // If al1 is better on frag1 and al2 on frag0: swap so al1 keeps frag0, al2 keeps frag1
            if (al1_frag0 - al2_frag0 < al1_frag1 - al2_frag1) {
                // al1 is relatively better on frag1, al2 on frag0 → swap
                swap(al1, al2);
            }

            // ex1: last exon of al1 on iFrag=0 (boundary before mate transition)
            ex1 = 0;
            for (uint i = 0; i < al1->nExons; i++) {
                if (al1->exons[i][EX_iFrag] == 0)
                    ex1 = i; // keeps updating to highest index on iFrag=0
            }
            // ex2: first exon of al2 on iFrag=1 (boundary after mate transition)
            ex2 = al2->nExons - 1; // default
            for (uint i = 0; i < al2->nExons; i++) {
                if (al2->exons[i][EX_iFrag] == 1) {
                    ex2 = i;
                    break;
                }
            }
        } else if (bestEx1 >= 0) {
            // Fall back to same-mate cross-chromosome overlap
            ex1 = bestEx1;
            ex2 = bestEx2;
        } else {
            // Fallback: positional selection (traditional chimeric junction between non-overlapping segments)
            ex1 = al1->Str==1 ? 0 : al1->nExons-1;
            ex2 = al2->Str==0 ? 0 : al2->nExons-1;
        }
    }
};

bool ChimericAlign::chimericCheck() {
    bool chimGood=true;

    chimGood = chimGood && al1->exons[ex1][EX_iFrag] <= al2->exons[ex2][EX_iFrag];//otherwise - strange configuration, both segments contain two mates
        //if ( trChim[0].exons[e0][EX_iFrag] > trChim[1].exons[e1][EX_iFrag] ) {//strange configuration, rare, similar to the next one
        //    chimN=0;//reject such chimeras
            //good test example:
            //CTTAGCTAGCAGCGTCTTCCCAGTGCCTGGAGGGCCAGTGAGAATGGCACCCTCTGGGATTTTTGCTCCTAGGTCT
            //TTGAGGTGAAGTTCAAAGATGTGGCTGGCTGTGAGGAGGCCGAGCTAGAGATCATGGAATTTGTGAATTTCTTGAA
        //} else

    //junction overhangs too short for chimerically spliced mates
    chimGood = chimGood && (al1->exons[ex1][EX_iFrag] < al2->exons[ex2][EX_iFrag] || (al1->exons[ex1][EX_L] >= pCh.junctionOverhangMin &&  al2->exons[ex2][EX_L] >= pCh.junctionOverhangMin) );

    return chimGood;
};
