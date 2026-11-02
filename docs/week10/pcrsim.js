/* PCR 시뮬레이션 엔진(브라우저 전용, 외부 전송 없음)
 *
 * 근거
 *  - 프라이머–주형 결합: nearest-neighbor 열역학(SantaLucia 1998), Mg2+ 보정(von Ahsen 2001).
 *    짝이 안 맞는 염기는 그 양옆 이웃 쌍을 약한 값으로 바꿔 계산(Allawi & SantaLucia 1997–98의 평균 크기).
 *  - 3′ 끝 불일치에 따른 연장 효율: Kwok et al. 1990 (Nucleic Acids Res), Stadhouders et al. 2010 (J Mol Diagn).
 *  - 비특이 결합 후보 찾기: 3′ 끝 짧은 조각이 맞는 자리를 먼저 찾고 전체 결합을 계산(MFEprimer, Qu et al. 2012의 방식).
 *  - 혼성화 속도: k_on ≈ 10^6 M⁻¹s⁻¹ (Wetmur 1991).
 *  - 중합효소 속도: Taq 약 60 nt/s (70 °C), 24 nt/s (55 °C) (Innis et al. 1988, PNAS)의 35%를 실효값으로 씀(72 °C에서 약 1.3 kb/분, "1 kb당 1분" 관행에 가까움).
 *  - 변성: 산물 Tm보다 충분히 뜨거워야 하고, 길수록 가닥이 다 떨어지는 데 시간이 더 걸린다고 봄.
 *  - 나머지 전사체(약 5×10^7 nt)·유전체(3×10^9 nt)의 무작위 결합: 3′ 끝 8–16 nt가 우연히 맞는 자리 수(4^-k)와
 *    그 결합 세기로 계산해, 결합 온도가 낮을 때 여러 크기 산물(번짐)이 생기게 함.
 *  - 사이클마다 증폭 효율이 프라이머 결합·연장·변성에 따라 정해지고, 산물이 쌓이면 효소·프라이머가 모자라 평탄해지는 모형
 *    (Mehra & Hu 2005; Gevertz et al. 2005, Biotechnol Bioeng).
 *  - 배경 서열: NCBI BLAST(refseq_rna, 사람)로 찾은 다른 유전자 전사체, 흔한 전사체, 유전체 구간(pcr-bg.js).
 */
(function (root) {
  var D = root.DNA, NN = D.NN, R = 1.987;
  var MIX = { mono: 50, mg: 3, dntp: 0.8, primer: 0.25 };   // SYBR qPCR 마스터믹스 수준
  var NAEQ = MIX.mono + 120 * Math.sqrt(Math.max(0, MIX.mg - MIX.dntp));
  var LNNA = Math.log(NAEQ / 1000), C = MIX.primer * 1e-6;
  var P0 = C * 20e-6 * 6.022e23;          // 20 µL 반응액의 프라이머 분자 수(약 3 × 10^12)
  var MMAX = 2.8e14;                      // 쌓일 수 있는 이중가닥 DNA 최대량(bp 합, 약 300 ng)
  var KON = 1e6;
  var COPIES = { cdna: 3000, gdna: 3000, contam: 30 };
  var MM = [0, -1.6];
  var COMP = { A: 'T', T: 'A', G: 'C', C: 'G' };
  var DEFAULT = { initT: 95, initS: 120, denT: 95, denS: 15, annT: 60, annS: 60, twoStep: true, extT: 72, extS: 30, cycles: 40 };

  function isB(b) { return b === 'A' || b === 'C' || b === 'G' || b === 'T'; }
  function initP(b) { return (b === 'G' || b === 'C') ? [0.1, -2.8] : [2.3, 4.1]; }

  // ── 결합 열역학 ─────────────────────────────────────────
  // p: 프라이머 5′→3′, s: 주형에서 프라이머와 같은 방향으로 읽은 자리(완전히 맞으면 s === p)
  function duplex(p, s) {
    var L = p.length, m = new Array(L), i, mis = 0;
    for (i = 0; i < L; i++) { m[i] = p[i] === s[i] && isB(p[i]); if (!m[i]) mis++; }
    var i1 = L - 1; while (i1 >= 0 && !m[i1]) i1--;
    if (i1 < 7) return null;
    // 3′ 쪽 끝(i1)에서 시작해 5′ 쪽으로 늘려 가며 가장 안정한 구간을 고른다.
    // 5′ 쪽에 틀린 염기가 몰려 있으면 그 부분은 붙지 않고 떠 있는 것으로 본다.
    var best = null, h = 0, sS = 0, Tref = 328.15;
    for (i = i1 - 1; i >= 0; i--) {
      var v = (m[i] && m[i + 1]) ? NN[p.substr(i, 2)] : MM; h += v[0]; sS += v[1];
      if (!m[i] || i1 - i < 7) continue;
      var a = initP(p[i]), b = initP(p[i1]);
      var dHi = h + a[0] + b[0], dSi = sS + a[1] + b[1] + 0.368 * (i1 - i) * LNNA, g = dHi - Tref * dSi / 1000;
      if (!best || g < best.g) best = { g: g, dH: dHi, dS: dSi };
    }
    if (!best) return null;
    var dH = best.dH, dS = best.dS;
    var mm3 = [];
    for (i = 1; i <= 5 && i <= L; i++) if (!m[L - i]) mm3.push(i);
    return { dH: dH, dS: dS, mis: mis, mm3: mm3, m: m };
  }
  // u: 머리핀으로 접히지 않은 프라이머 비율. 접힌 프라이머만큼 결합에 쓸 수 있는 농도가 줄어든다
  function theta(dH, dS, T, u) {
    var TK = T + 273.15, K = Math.exp(-(dH * 1000 - TK * dS) / (R * TK)), Ce = C * (u == null ? 1 : u);
    return K * Ce / (1 + K * Ce);
  }
  // 반응 조건(프라이머가 주형보다 훨씬 많음)에서 절반이 붙는 온도
  function tmSite(dH, dS) { return dH * 1000 / (dS + R * Math.log(C)) - 273.15; }

  // 3′ 끝에서 k번째 염기가 안 맞을 때 연장 효율(프라이머:주형 염기 종류에 따라)
  var EXT3 = { 1: [0.03, 0.005, 0.3], 2: [0.15, 0.05, 0.5], 3: [0.4, 0.25, 0.7], 4: [0.65, 0.5, 0.85], 5: [0.85, 0.75, 0.95] };
  function ext3(p, s, mm3) {
    var f = 1;
    mm3.forEach(function (k) {
      var i = p.length - k, pb = p[i], tb = COMP[s[i]] || 'N';
      var strong = (pb === 'A' && tb === 'G') || (pb === 'G' && tb === 'A') || (pb === 'C' && tb === 'C');
      f *= strong ? EXT3[k][1] : pb === 'T' ? EXT3[k][2] : EXT3[k][0];
    });
    return f;
  }

  // 머리핀: 접힌 비율
  function hairpinDG(s) {
    var h = D.hairpin(s);
    if (!h.stem || h.stem < 3) return null;
    var dH = 0, dS = 0;
    for (var x = 0; x < h.stem - 1; x++) { var v = NN[s.substr(h.i + x, 2)]; if (!v) return null; dH += v[0]; dS += v[1]; }
    var loopG = h.loop <= 4 ? 3.5 : h.loop <= 6 ? 3.3 : 3.3 + 1.75 * 0.616 * Math.log(h.loop / 6);
    dS += -loopG * 1000 / 310.15 + 0.368 * (h.stem - 1) * LNNA;
    return { dH: dH, dS: dS, stem: h.stem, loop: h.loop, end3: h.end3 };
  }
  function dG(dH, dS, T) { return dH - (T + 273.15) * dS / 1000; }
  function unfolded(s, T) {
    var hp = hairpinDG(s); if (!hp) return 1;
    var TK = T + 273.15, K = Math.exp(-(hp.dH * 1000 - TK * hp.dS) / (R * TK));
    return 1 / (1 + K);
  }

  // 프라이머 이합체: a의 3′ 끝 또는 b의 3′ 끝이 걸린 가장 긴 짝
  function dimerInfo(a, b, T) {
    var d = D.dimer(a, b).best3;
    if (!d.run || d.run < 2) return null;
    var seg = a.substr(d.aStart, d.run), dH = 0, dS = 0;
    for (var x = 0; x < seg.length - 1; x++) { var v = NN[seg.substr(x, 2)]; dH += v[0]; dS += v[1]; }
    var i1 = initP(seg[0]), i2 = initP(seg[seg.length - 1]);
    dH += i1[0] + i2[0]; dS += i1[1] + i2[1] + 0.368 * (seg.length - 1) * LNNA;
    var out = { run: d.run, dH: dH, dS: dS, dG: dG(dH, dS, T), theta: theta(dH, dS, T), size: a.length + b.length - d.run, info: d };
    // 틀린 염기가 섞인 겹침까지: a의 3′ 끝이 b 위 어디에 얹혀도 되고, b가 a의 3′ 끝 너머로 이어져야 연장된다
    var brc = D.rc(b), La = a.length, Lb = b.length, bestP = out.theta;
    for (var q = 0; q < Lb - 1; q++) {
      var s = '';
      for (var i = 0; i < La; i++) { var j = q - (La - 1 - i); s += j >= 0 ? brc[j] : 'N'; }
      var dd = duplex(a, s);
      if (!dd) continue;
      var th = theta(dd.dH, dd.dS, T) * ext3(a, s, dd.mm3);
      if (th > bestP) { bestP = th; out.size = La + Lb - 1 - q; }
    }
    out.theta = bestP;
    return out;
  }

  // 중합효소 실효 속도(nt/s)
  var RATE = [[30, 0.5], [37, 1.5], [45, 6], [55, 24], [60, 35], [65, 48], [70, 60], [75, 65], [80, 55], [85, 30], [90, 8], [95, 0]];
  function rate(T) {
    if (T <= RATE[0][0]) return RATE[0][1] * RATE_EFF;
    for (var i = 1; i < RATE.length; i++) if (T <= RATE[i][0]) {
      var a = RATE[i - 1], b = RATE[i]; return (a[1] + (b[1] - a[1]) * (T - a[0]) / (b[0] - a[0])) * RATE_EFF;
    }
    return 0;
  }
  var RATE_EFF = 0.35;   // 72 °C에서 약 21 nt/s (1.3 kb/분). 흔히 쓰는 "1 kb당 1분" 기준에 가깝게

  // 변성: 산물 Tm보다 충분히 뜨거워야 하고, 길수록 다 떨어지는 데 시간이 더 걸린다
  function denFrac(L, tm, prog) {
    var tau = 0.4 + L / 250;
    return (1 / (1 + Math.exp((tm + 1 - prog.denT) / 1.5))) * (1 - Math.exp(-(prog.denS + 0.3) / tau));
  }
  function extFrac(L, dist) { return 1 / (1 + Math.pow(L / Math.max(1, dist), 6)); }

  // ── 나머지 전사체·유전체에 대한 무작위 결합(통계 근사) ──
  // 실제 서열을 다 넣을 수 없으므로, 크기 B인 무작위 서열에서 프라이머 3′ 끝 k nt가 우연히 맞는 자리 수
  // n_k = 2B · (3/4) · 4^-k 와 그 k-mer의 결합 비율 θ_k를 곱해 더한다. 결합 온도가 낮을수록 짧은 일치도 붙는다.
  // cDNA: 서로 다른 전사체 서열 약 5×10^7 nt, 반응 하나에 전사체 약 2×10^8개(평균 2 kb) → 같은 자리 평균 8,000 copies
  var BG = { cdna: { B: 5e7, c: 8000 }, gdna: { B: 3e9 } };
  // 3′ 끝 k nt가 맞고 그 앞은 무작위(4분의 1 확률로 맞음)인 자리를 뽑아 결합 비율의 평균을 낸다.
  // 결합 온도가 낮으면 틀린 염기가 몇 개 섞인 자리도 붙기 때문에 k가 짧은 자리가 크게 늘어난다.
  var BASES = 'ACGT';
  // 서열 구성: 사람 유전체는 GC 41%, CG 이웃 쌍이 기대값의 약 1/4 (CpG 결핍). 전사체는 GC 약 50%, CG 약 1/2로 둔다.
  var COMPO = { cdna: { gc: 0.50, cpg: 0.5 }, gdna: { gc: 0.41, cpg: 0.24 } };
  function kmerFreq(s, c) {
    var f = 1;
    for (var i = 0; i < s.length; i++) f *= (s[i] === 'G' || s[i] === 'C') ? c.gc / 2 : (1 - c.gc) / 2;
    var cg = (s.match(/CG/g) || []).length;
    return f * Math.pow(c.cpg, cg);
  }
  function primeP(p, T, unf, kin, B, compo) {
    var P = 0, L = p.length, h = 0;
    for (var i = 0; i < L; i++) h = (h * 31 + p.charCodeAt(i)) >>> 0;
    var rr = rng(h || 1), NS = 80;
    for (var k = 6; k <= Math.min(16, L - 1); k++) {
      var nk = 2 * B * 0.75 * kmerFreq(p.substr(L - k), compo), sum = 0;
      for (var n = 0; n < NS; n++) {
        var s = '';
        for (var j = 0; j < L; j++) {
          if (j >= L - k) s += p[j];
          else if (j === L - k - 1) { var b; do { b = BASES[(rr() * 4) | 0]; } while (b === p[j]); s += b; }
          else s += rr() < 0.25 ? p[j] : (function () { var b2; do { b2 = BASES[(rr() * 4) | 0]; } while (b2 === p[j]); return b2; })();
        }
        var d = duplex(p, s);
        if (d) sum += theta(d.dH, d.dS, T, unf);
      }
      P += nk * sum / NS;
    }
    return P * kin;
  }
  function meltTm(gc, n) { return 81.5 + 0.41 * gc - 675 / n + 16.6 * Math.log(0.15) / Math.LN10; }

  // ── 난수(실행마다 조금씩 다른 결과) ─────────────────────
  function rng(seed) {
    var a = seed >>> 0;
    return function () { a = (a + 0x6D2B79F5) >>> 0; var t = a; t = Math.imul(t ^ (t >>> 15), t | 1); t ^= t + Math.imul(t ^ (t >>> 7), t | 61); return ((t ^ (t >>> 14)) >>> 0) / 4294967296; };
  }
  function gauss(r) { var u = 1 - r(), v = r(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v); }
  function poisson(lam, r) {
    if (lam <= 0) return 0;
    if (lam > 40) return Math.max(0, Math.round(lam + Math.sqrt(lam) * gauss(r)));
    var L = Math.exp(-lam), k = 0, p = 1;
    do { k++; p *= r(); } while (p > L);
    return k - 1;
  }

  // ── 결합 자리 찾기 ───────────────────────────────────────
  // 3′ 끝 6 nt, 또는 마지막 염기를 뺀 6 nt가 그대로 맞는 곳을 후보로 잡고 전체 결합을 계산
  function scan(strand, p, out, dir, tlen, T, kin, unf) {
    var L = p.length, seeds = [[p.substr(L - 6), L - 6], [p.substr(L - 7, 6), L - 7]], seen = {};
    seeds.forEach(function (sd) {
      for (var i = strand.indexOf(sd[0]); i >= 0; i = strand.indexOf(sd[0], i + 1)) {
        var st = i - sd[1];
        if (st < 0 || st + L > strand.length || seen[st]) continue;
        seen[st] = 1;
        var s = strand.substr(st, L), d = duplex(p, s);
        if (!d || d.mis > 7) continue;
        var th = theta(d.dH, d.dS, T, unf), e3 = ext3(p, s, d.mm3), phi = th * e3 * kin;
        if (phi < 1e-9 && d.mis > 0) continue;
        var top = dir === '+' ? st : tlen - st - L;
        out.push({ start: top, end: top + L - 1, dir: dir, mis: d.mis, mm3: d.mm3, theta: th, e3: e3, phi: phi, tm: tmSite(d.dH, d.dS), site: s, m: d.m });
      }
    });
  }
  function sitesOn(t, trc, p, T, kin, unf) {
    var out = [];
    scan(t, p, out, '+', t.length, T, kin, unf);
    scan(trc, p, out, '-', t.length, T, kin, unf);
    return out;
  }

  function bgFor(gene) { return (root.PCR_BG && gene && root.PCR_BG.genes[gene]) || null; }
  function detectGene(raw, M) {
    var bg = root.PCR_BG; if (!bg) return null;
    var head = String(raw || '').split(/\r?\n/)[0], acc = (head.match(/[NX][MR]_\d+/) || [])[0];
    for (var g in bg.genes) {
      var x = bg.genes[g];
      if (acc && x.nm.split('.')[0] === acc) return g;
      if (M && M.substr(0, 40) === x.sig) return g;
    }
    return null;
  }

  // ── 시뮬레이션 ───────────────────────────────────────────
  // opts: { M, G, gene, contam, prog, pairs:[{name,F,R}], seed }
  function run(opts) {
    var prog = Object.assign({}, DEFAULT, opts.prog || {});
    var r = rng(opts.seed == null ? (Date.now() & 0x7fffffff) : opts.seed);
    var bg = bgFor(opts.gene), common = (root.PCR_BG && root.PCR_BG.common) || [];
    var T = prog.annT, kin = 1 - Math.exp(-(prog.annS + 1) * KON * C);
    // 한 사이클에 중합효소가 갈 수 있는 거리: 결합 단계(붙은 뒤 남은 시간의 일부) + 연장 단계 + 온도 올라가는 사이
    var extS = prog.twoStep ? 0 : prog.extS;
    var dist = rate(T) * prog.annS * (prog.twoStep || !extS ? 0.8 : 0.3) + rate(prog.extT) * extS + (extS ? rate((T + prog.extT) / 2) * 0.5 : 0) + rate((T + prog.denT) / 2) * 0.3;

    function lanesFor() {
      var cd = [{ name: opts.geneLabel || '내 mRNA', kind: 'target', seq: opts.M, copies: COPIES.cdna }];
      common.forEach(function (x) { cd.push({ name: x.name, acc: x.acc, kind: 'common', seq: x.seq, copies: x.copies * COPIES.cdna }); });
      if (bg) bg.cdna.forEach(function (x) { cd.push({ name: x.name, acc: x.acc, kind: x.kind, seq: x.seq, copies: x.copies * COPIES.cdna, ident: x.ident }); });
      var gen = [];
      if (opts.G) {
        gen.push({ name: (opts.geneLabel || '내 유전자') + ' 유전체', kind: 'gdna', seq: opts.G, copies: 1 });
        if (bg) bg.gdna.forEach(function (x) {
          var seq = x.seq;
          if (!seq && x.ref) bg.cdna.forEach(function (y) { if (y.name === x.ref) seq = y.seq; });
          if (seq) gen.push({ name: x.name + (x.kind === 'repeat' ? '' : ' 유전체'), acc: x.acc, kind: x.kind === 'repeat' ? 'repeat' : 'gdna-other', seq: seq, copies: x.copies });
        });
      }
      var L = [];
      var cdna = cd.slice();
      if (opts.contam) gen.forEach(function (x) { cdna.push(Object.assign({}, x, { copies: x.copies * COPIES.contam, contam: true })); });
      L.push({ name: 'cDNA', templates: cdna });
      if (opts.G) L.push({ name: 'gDNA', templates: gen.map(function (x) { return Object.assign({}, x, { copies: x.copies * COPIES.gdna }); }) });
      L.push({ name: 'NTC', templates: [] });
      return L;
    }

    var out = { prog: prog, dist: dist, kin: kin, pairs: [] };
    (opts.pairs || []).forEach(function (pair) {
      var prs = [{ name: 'F', seq: pair.F }, { name: 'R', seq: pair.R }];
      prs.forEach(function (p) {
        p.unf = unfolded(p.seq, T);
        var d = duplex(p.seq, p.seq);
        p.dH = d.dH; p.dS = d.dS; p.theta = theta(d.dH, d.dS, T, p.unf); p.tmMix = tmSite(d.dH, d.dS);
        p.phi = p.theta * kin;
        var pb = D.tmNN(p.seq, D.DESIGN); p.tm = pb ? pb.tm : NaN;
        p.Pc = primeP(p.seq, T, p.unf, kin, BG.cdna.B, COMPO.cdna);
        p.Pg = primeP(p.seq, T, p.unf, kin, BG.gdna.B, COMPO.gdna);
      });
      var PcT = prs[0].Pc + prs[1].Pc, PgT = prs[0].Pg + prs[1].Pg;
      // 프라이머 이합체 후보
      var dims = [];
      [[0, 1], [1, 0], [0, 0], [1, 1]].forEach(function (ij) {
        var a = prs[ij[0]], b = prs[ij[1]], di = dimerInfo(a.seq, b.seq, T);
        if (!di) return;
        var lam = P0 * di.theta * 1e-6;
        dims.push({ X: a, Y: b, label: a.name + '+' + b.name, run: di.run, dGT: di.dG, lam: lam, size: Math.max(di.size, 25), gc: D.gcPercent(a.seq + b.seq), info: di.info });
      });

      var lanes = lanesFor().map(function (lane) {
        var products = [];
        lane.templates.forEach(function (tp) {
          var t = tp.seq, trc = D.rc(t), S = {};
          prs.forEach(function (p) { S[p.name] = sitesOn(t, trc, p.seq, T, kin, p.unf); });
          tp.sites = S;
          prs.forEach(function (X) {
            S[X.name].filter(function (a) { return a.dir === '+'; }).forEach(function (a) {
              prs.forEach(function (Y) {
                S[Y.name].filter(function (b) { return b.dir === '-'; }).forEach(function (b) {
                  var size = b.end - a.start + 1;
                  if (size <= Math.max(X.seq.length, Y.seq.length) || size > 6000) return;
                  var seedPhi = a.phi * b.phi;
                  if (seedPhi * tp.copies < 1e-7) return;
                  var gc = D.gcPercent(t.substring(a.start, b.end + 1)), tm = meltTm(gc, size);
                  products.push({ X: X, Y: Y, a: a, b: b, size: size, gc: gc, tm: tm, tp: tp, T0: tp.copies, seedPhi: seedPhi, N: 0 });
                });
              });
            });
          });
        });
        // 무작위 결합에서 생기는 여러 크기의 산물(번짐): 사이클마다 새로 생기는 수
        var reach = Math.min(dist, 3000), sm = [];
        if (reach > 80) {
          var rate0 = 0;
          if (lane.name === 'cDNA') {
            rate0 = BG.cdna.c * PcT * PcT * (reach - 60) / (4 * BG.cdna.B);
            if (opts.contam && opts.G) rate0 += COPIES.contam * PgT * PgT * (reach - 60) / (4 * BG.gdna.B);
          } else if (lane.name === 'gDNA') rate0 = COPIES.gdna * PgT * PgT * (reach - 60) / (4 * BG.gdna.B);
          if (rate0 > 1e-6) {
            var nb = 10, lo = Math.log(60), hi = Math.log(reach);
            for (var q = 0; q < nb; q++) {
              var a1 = Math.exp(lo + (hi - lo) * q / nb), a2 = Math.exp(lo + (hi - lo) * (q + 1) / nb), L = Math.round(Math.sqrt(a1 * a2));
              sm.push({ size: L, w: (a2 - a1) / (reach - 60), rate: rate0 * (a2 - a1) / (reach - 60), tm: meltTm(48, L), N: 0 });
            }
          }
        }
        return { name: lane.name, templates: lane.templates, products: products, smearP: sm };
      });

      lanes.forEach(function (lane) {
        var noise = Math.exp(0.1 * gauss(r));
        var pool = { F: P0, R: P0 }, Mtot = 0, curve = [], ps = lane.products, sm = lane.smearP;
        ps.forEach(function (p) {
          var c = extFrac(p.size, dist), den = denFrac(p.size, p.tm, prog);
          var Emax = 0.97 * (p.gc > 65 ? 0.9 : 1);
          p.c = c; p.den = den;
          p.E = Emax * Math.sqrt(p.X.phi * p.Y.phi) * c * den;
          p.seed = p.T0 * noise * p.seedPhi * c;
          p.genomic = p.tp.kind === 'gdna' || p.tp.kind === 'gdna-other' || p.tp.kind === 'repeat';
        });
        var ePerf = Math.sqrt(prs[0].phi * prs[1].phi);
        sm.forEach(function (x) { x.E = 0.95 * ePerf * extFrac(x.size, dist) * denFrac(x.size, x.tm, prog); });
        var ds = dims.map(function (d) { var tm = meltTm(d.gc, d.size); return { d: d, N: 0, E: 0.9 * Math.sqrt(d.X.phi * d.Y.phi) * extFrac(d.size, dist) * denFrac(d.size, tm, prog), size: d.size, tm: tm, X: d.X, Y: d.Y }; });
        for (var k = 1; k <= prog.cycles; k++) {
          var A = 1 - Math.exp(-(prog.initS + (k - 1) * prog.denS) / 40);
          var s = Math.max(0, 1 - Mtot / MMAX);
          var gDen = 1 - Math.exp(-(k === 1 ? prog.initS : prog.denS) / 10);
          var gains = ps.map(function (p) {
            var pf = Math.min(pool[p.X.name], pool[p.Y.name]) / P0;
            return (p.N * p.E + p.seed * (p.genomic ? gDen : 1)) * A * s * pf;
          });
          var dg = ds.map(function (x) {
            var pf = Math.min(pool[x.X.name], pool[x.Y.name]) / P0;
            return x.N * x.E * A * s * pf + poisson(x.d.lam * A * s * pf, r);
          });
          var pfs = Math.min(pool.F, pool.R) / P0;
          var sg = sm.map(function (x) { return (x.N * x.E + x.rate * noise * (k === 1 ? gDen : 1)) * A * s * pfs; });
          ps.forEach(function (p, i) { p.N += gains[i]; pool[p.X.name] -= gains[i]; pool[p.Y.name] -= gains[i]; });
          ds.forEach(function (x, i) { x.N += dg[i]; pool[x.X.name] -= dg[i]; pool[x.Y.name] -= dg[i]; });
          sm.forEach(function (x, i) { x.N += sg[i]; pool.F -= sg[i]; pool.R -= sg[i]; });
          pool.F = Math.max(0, pool.F); pool.R = Math.max(0, pool.R);
          Mtot = 0;
          ps.forEach(function (p) { Mtot += p.N * p.size; });
          ds.forEach(function (x) { Mtot += x.N * x.size; });
          sm.forEach(function (x) { Mtot += x.N * x.size; });
          curve.push(Mtot);
        }
        // 형광(ΔRn): 이중가닥 DNA 양에 비례 + 바탕 잡음
        var base = 0.012 * (1 + 0.3 * r());
        lane.fluor = curve.map(function (m) { return m / 1e12 + base * (1 + 0.15 * gauss(r)); });
        lane.ct = null;
        for (var j = 0; j < lane.fluor.length; j++) if (lane.fluor[j] >= THRESH) {
          var f0 = j ? lane.fluor[j - 1] : base, f1 = lane.fluor[j];
          lane.ct = j + (Math.log(THRESH) - Math.log(Math.max(f0, 1e-6))) / (Math.log(f1) - Math.log(Math.max(f0, 1e-6))); // 사이클은 1부터
          break;
        }
        lane.dimers = ds;
        lane.total = Mtot;
        // 젤에 올리는 양: 20 µL 중 5 µL
        lane.bands = [];
        ps.forEach(function (p) { if (p.N > 0) lane.bands.push({ size: p.size, ng: p.N * p.size * 650 / 6.022e23 * 1e9 * 0.25, tm: p.tm, p: p }); });
        ds.forEach(function (x) { if (x.N > 0) lane.bands.push({ size: x.size, ng: x.N * x.size * 650 / 6.022e23 * 1e9 * 0.25, tm: x.tm, dimer: x }); });
        lane.bands = mergeBands(lane.bands);
        lane.smear = sm.map(function (x) { return { size: x.size, ng: x.N * x.size * 650 / 6.022e23 * 1e9 * 0.25, tm: x.tm }; });
        lane.smearNg = lane.smear.reduce(function (a, b) { return a + b.ng; }, 0);
      });
      out.pairs.push({ name: pair.name, F: pair.F, R: pair.R, primers: prs, dimers: dims, lanes: lanes });
    });
    return out;
  }
  var THRESH = 1.5;

  // 크기가 3% 안으로 같은 산물은 젤에서 한 밴드
  function mergeBands(bs) {
    bs.sort(function (a, b) { return b.size - a.size; });
    var out = [];
    bs.forEach(function (b) {
      var o = out.find(function (x) { return Math.abs(x.size - b.size) / x.size < 0.03 && !!x.dimer === !!b.dimer; });
      if (o) { o.parts.push(b); if (b.ng > o.main.ng) o.main = b; o.ng += b.ng; }
      else out.push({ size: b.size, ng: b.ng, tm: b.tm, dimer: !!b.dimer, parts: [b], main: b });
    });
    out.forEach(function (o) { o.size = o.main.size; o.tm = o.main.tm; });
    return out;
  }

  root.PCRSIM = {
    run: run, DEFAULT: DEFAULT, MIX: MIX, COPIES: COPIES, THRESH: THRESH, detectGene: detectGene,
    duplex: duplex, theta: theta, tmSite: tmSite, unfolded: unfolded, hairpinDG: hairpinDG, dimerInfo: dimerInfo, rate: rate, meltTm: meltTm, dG: dG,
    VISIBLE_NG: 0.5
  };
})(window);
