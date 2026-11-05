/* r1_geom_scan.c -- classify every Flex v2 R1 read the way piscem-rs 0.11.0 parses it,
 * and record why a read is dropped before mapping.
 *
 * The geometry comes from the command line and follows
 * GEOM_V2 = 1{b[CB]u[UMI]x[LO-HI]hamming(f[ANCHOR],D)s[SBC]x:}2{r:}:
 *
 *   window   anchor start offsets CB+UMI+LO ..= CB+UMI+HI
 *   anchor   seq_geom_parser 1.1.0 find_anchor(): scan the window left to right, keep the
 *            first position of minimum Hamming distance, stop early at distance 0; accept if
 *            the distance is <= D.  Bytes are compared raw, so an N is a mismatch.
 *   fields   piscem-rs 0.11.0 processors.rs (multi-barcode path): the UMI must be all ACGT;
 *            each barcode (b0 = sample, then b1 = cell) may hold at most one N (recovered to
 *            A) and must otherwise be ACGT.  A field running past the read end is missing.
 *
 * A dropped read is counted in num_reads but in neither num_mapped, num_poisoned nor
 * unmapped_bc_count.bin, so  sum(dropped) == num_reads - num_mapped - num_poisoned - ledger.
 *
 * Drop precedence (one class per read): anchor_absent_in_window, umi_invalid,
 * sample_bc_invalid, cell_bc_invalid.  For anchor-absent reads the whole read is searched for
 * the anchor (first exact occurrence; failing that, first occurrence within distance D) and
 * the read is sub-classified by where that anchor sits relative to the window.
 *
 * usage: r1_geom_scan --cb-len 16 --umi-len 12 --win-lo 0 --win-hi 4 --anchor TTGCTAGGACCG
 *                     --max-dist 1 --sbc-len 10 --sbc-list sample_bc_list.tsv  < R1.fastq
 * Reads uncompressed FASTQ on stdin; prints a JSON summary on stdout.
 * cc -O2 -o r1_geom_scan r1_geom_scan.c
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXOFF 256
#define MAXSBC 64

static int CB_LEN, UMI_LEN, WIN_LO, WIN_HI, ALEN, MAXD, SBC_LEN;
static char ANCHOR[64];
static char SBC[MAXSBC][64];
static int NSBC;

enum { C_PARSE_OK = 0, C_ANCHOR_MISS, C_UMI_BAD, C_SBC_INVALID, C_CBC_INVALID, N_CLASS };
static const char *CLASS_NAME[N_CLASS] = {
    "parse_ok", "anchor_absent_in_window", "umi_invalid", "sample_bc_invalid", "cell_bc_invalid"};

enum { A_BEFORE = 0, A_AFTER, A_NONE, N_ASUB };
static const char *ASUB_NAME[N_ASUB] = {
    "anchor_before_window", "anchor_after_window", "no_anchor_within_max_dist"};

static inline int is_n(char c) { return c == 'N' || c == 'n'; }
static inline int is_acgt(char c) {
    return c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'a' || c == 'c' || c == 'g' || c == 't';
}
static inline int all_acgt(const char *p, int len) {
    for (int i = 0; i < len; i++) if (!is_acgt(p[i])) return 0;
    return 1;
}
/* piscem barcode validity: <= 1 N, every other base ACGT */
static inline int bc_valid(const char *p, int len) {
    int nn = 0;
    for (int i = 0; i < len; i++) {
        if (is_n(p[i])) { if (++nn > 1) return 0; }
        else if (!is_acgt(p[i])) return 0;
    }
    return 1;
}
static inline int ham_at_most(const char *a, const char *b, int len, int limit) {
    int d = 0;
    for (int i = 0; i < len; i++) if (a[i] != b[i] && ++d > limit) return d;
    return d;
}
/* seq_geom_parser find_anchor over [lo, hi]; returns position, sets *dist */
static int find_anchor(const char *s, int L, int lo, int hi, int *dist) {
    int best = MAXD + 1, bestp = -1;
    for (int p = lo; p <= hi; p++) {
        if (p + ALEN > L) continue;
        int d = ham_at_most(s + p, ANCHOR, ALEN, MAXD);
        if (d < best) { best = d; bestp = p; if (d == 0) break; }
    }
    *dist = best;
    return best <= MAXD ? bestp : -1;
}
static int sbc_index(const char *p, int avail) {
    if (avail < SBC_LEN) return -1;
    for (int i = 0; i < NSBC; i++) if (memcmp(p, SBC[i], SBC_LEN) == 0) return i;
    return -1;
}

static void die(const char *m) { fprintf(stderr, "r1_geom_scan: %s\n", m); exit(2); }

int main(int argc, char **argv) {
    const char *sbc_file = NULL;
    CB_LEN = UMI_LEN = WIN_LO = WIN_HI = MAXD = SBC_LEN = -1;
    for (int i = 1; i < argc; i++) {
        if (i + 1 >= argc) die("option without value");
        const char *k = argv[i], *v = argv[++i];
        if (!strcmp(k, "--cb-len")) CB_LEN = atoi(v);
        else if (!strcmp(k, "--umi-len")) UMI_LEN = atoi(v);
        else if (!strcmp(k, "--win-lo")) WIN_LO = atoi(v);
        else if (!strcmp(k, "--win-hi")) WIN_HI = atoi(v);
        else if (!strcmp(k, "--max-dist")) MAXD = atoi(v);
        else if (!strcmp(k, "--sbc-len")) SBC_LEN = atoi(v);
        else if (!strcmp(k, "--anchor")) { if (strlen(v) >= sizeof ANCHOR) die("anchor too long"); strcpy(ANCHOR, v); }
        else if (!strcmp(k, "--sbc-list")) sbc_file = v;
        else die("unknown option");
    }
    if (CB_LEN < 1 || UMI_LEN < 1 || WIN_LO < 0 || WIN_HI < WIN_LO || MAXD < 0 || SBC_LEN < 1 ||
        !ANCHOR[0] || !sbc_file) die("missing or invalid option");
    ALEN = (int)strlen(ANCHOR);
    /* window is given as spacer lengths; convert to absolute anchor offsets */
    WIN_LO += CB_LEN + UMI_LEN;
    WIN_HI += CB_LEN + UMI_LEN;

    FILE *sf = fopen(sbc_file, "r");
    if (!sf) die("cannot open --sbc-list");
    char line[1024];
    while (fgets(line, sizeof line, sf)) {
        char *t = strtok(line, "\t\r\n");
        if (!t || !*t) continue;
        if ((int)strlen(t) != SBC_LEN) die("sample barcode length != --sbc-len");
        if (NSBC == MAXSBC) die("too many sample barcodes");
        strcpy(SBC[NSBC++], t);
    }
    fclose(sf);
    if (!NSBC) die("empty --sbc-list");

    static char buf[1 << 20];
    uint64_t nreads = 0, cls[N_CLASS] = {0}, asub[N_ASUB] = {0};
    /* retained reads by (anchor offset, distance) */
    static uint64_t ok_off_dist[MAXOFF][8];
    /* anchor-absent reads by offset of the anchor found outside the window, and its distance */
    static uint64_t out_off_dist[MAXOFF][8];
    uint64_t ok_sbc_listed = 0, ok_sbc_unlisted = 0;
    uint64_t ok_per_sample[MAXSBC] = {0};

    int ln = 0;
    while (fgets(buf, sizeof buf, stdin)) {
        if ((ln++ & 3) != 1) continue;
        int L = (int)strlen(buf);
        while (L > 0 && (buf[L - 1] == '\n' || buf[L - 1] == '\r')) L--;
        buf[L] = '\0';
        const char *s = buf;
        nreads++;

        int d = 0;
        int pos = find_anchor(s, L, WIN_LO, WIN_HI, &d);
        int klass;
        if (pos < 0) klass = C_ANCHOR_MISS;
        else if (CB_LEN + UMI_LEN > L || !all_acgt(s + CB_LEN, UMI_LEN)) klass = C_UMI_BAD;
        else if (pos + ALEN + SBC_LEN > L || !bc_valid(s + pos + ALEN, SBC_LEN)) klass = C_SBC_INVALID;
        else if (!bc_valid(s, CB_LEN)) klass = C_CBC_INVALID;
        else klass = C_PARSE_OK;
        cls[klass]++;

        if (klass == C_PARSE_OK) {
            ok_off_dist[pos < MAXOFF ? pos : MAXOFF - 1][d < 8 ? d : 7]++;
            int si = sbc_index(s + pos + ALEN, L - (pos + ALEN));
            if (si >= 0) { ok_sbc_listed++; ok_per_sample[si]++; }
            else ok_sbc_unlisted++;
        } else if (klass == C_ANCHOR_MISS) {
            int dd = 0;
            int p = -1;
            /* first exact occurrence anywhere; else first within MAXD anywhere */
            for (int q = 0; q + ALEN <= L; q++)
                if (memcmp(s + q, ANCHOR, ALEN) == 0) { p = q; dd = 0; break; }
            if (p < 0 && MAXD > 0) {
                for (int q = 0; q + ALEN <= L; q++) {
                    int e = ham_at_most(s + q, ANCHOR, ALEN, MAXD);
                    if (e <= MAXD) { p = q; dd = e; break; }
                }
            }
            if (p < 0) asub[A_NONE]++;
            else {
                asub[p < WIN_LO ? A_BEFORE : A_AFTER]++;
                out_off_dist[p < MAXOFF ? p : MAXOFF - 1][dd < 8 ? dd : 7]++;
            }
        }
    }

    printf("{\n  \"num_reads\": %llu,\n", (unsigned long long)nreads);
    printf("  \"geometry\": {\"cb_len\": %d, \"umi_len\": %d, \"anchor\": \"%s\", \"max_dist\": %d, "
           "\"window_abs_lo\": %d, \"window_abs_hi\": %d, \"sbc_len\": %d, \"n_sample_bcs\": %d},\n",
           CB_LEN, UMI_LEN, ANCHOR, MAXD, WIN_LO, WIN_HI, SBC_LEN, NSBC);
    printf("  \"classes\": {");
    for (int i = 0; i < N_CLASS; i++)
        printf("%s\"%s\": %llu", i ? ", " : "", CLASS_NAME[i], (unsigned long long)cls[i]);
    printf("},\n  \"anchor_absent_subclasses\": {");
    for (int i = 0; i < N_ASUB; i++)
        printf("%s\"%s\": %llu", i ? ", " : "", ASUB_NAME[i], (unsigned long long)asub[i]);
    printf("},\n");
    const char *hn[2] = {"parse_ok_by_offset_dist", "anchor_absent_by_found_offset_dist"};
    for (int h = 0; h < 2; h++) {
        uint64_t (*H)[8] = h ? out_off_dist : ok_off_dist;
        printf("  \"%s\": {", hn[h]);
        int first = 1;
        for (int o = 0; o < MAXOFF; o++)
            for (int k = 0; k < 8; k++)
                if (H[o][k]) {
                    printf("%s\"%d:%d\": %llu", first ? "" : ", ", o, k, (unsigned long long)H[o][k]);
                    first = 0;
                }
        printf("},\n");
    }
    printf("  \"parse_ok_sample_bc_exact_listed\": %llu,\n  \"parse_ok_sample_bc_not_listed\": %llu,\n",
           (unsigned long long)ok_sbc_listed, (unsigned long long)ok_sbc_unlisted);
    printf("  \"parse_ok_per_sample_exact\": [");
    for (int i = 0; i < NSBC; i++) printf("%s%llu", i ? ", " : "", (unsigned long long)ok_per_sample[i]);
    printf("]\n}\n");
    return 0;
}
