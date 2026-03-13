#include "ChimericAlign.h"
#include <iostream>

ChimericAlign::ChimericAlign(ChimericSegment &seg1in, ChimericSegment &seg2in, int chimScoreIn, Genome &genomeIn, ReadAlign *RAin)
                              : seg1(seg1in), seg2(seg2in),chimScore(chimScoreIn), RA(RAin), P(seg1in.P), pCh(P.pCh), mapGen(genomeIn) {
    stitchingDone=false;

    al1=&seg1.align;
    al2=&seg2.align;

    if (al1->roStart > al2->roStart)
        swap (al1,al2);

    ex1 = al1->Str==1 ? 0 : al1->nExons-1;

    // Compute ro range of al1's junction exon
    uint al1_juncRoS, al1_juncRoE;
    if (al1->Str==0) {
        al1_juncRoS = al1->exons[ex1][EX_R];
        al1_juncRoE = al1->exons[ex1][EX_R] + al1->exons[ex1][EX_L] - 1;
    } else {
        al1_juncRoS = al1->Lread - al1->exons[ex1][EX_R] - al1->exons[ex1][EX_L];
        al1_juncRoE = al1_juncRoS + al1->exons[ex1][EX_L] - 1;
    }

    // Find the exon of al2 that overlaps with al1's junction exon in ro-space
    ex2 = al2->Str==0 ? 0 : al2->nExons-1; // default (junction-facing end)
    for (uint iex=0; iex<al2->nExons; iex++) {
        uint roS, roE;
        if (al2->Str==0) {
            roS = al2->exons[iex][EX_R];
            roE = al2->exons[iex][EX_R] + al2->exons[iex][EX_L] - 1;
        } else {
            roS = al2->Lread - al2->exons[iex][EX_R] - al2->exons[iex][EX_L];
            roE = roS + al2->exons[iex][EX_L] - 1;
        }
        if (roS <= al1_juncRoE && roE >= al1_juncRoS) {
            ex2 = iex;
            break;
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
