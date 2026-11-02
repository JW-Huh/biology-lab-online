/* 프라이머 설계 실습 공통 계산 함수(브라우저 전용, 외부 전송 없음) */
(function (root) {
  var COMP = { A: 'T', T: 'A', G: 'C', C: 'G', N: 'N', R: 'Y', Y: 'R', S: 'S', W: 'W', K: 'M', M: 'K', B: 'V', V: 'B', D: 'H', H: 'D' };

  // FASTA 첫 줄, GenBank의 ORIGIN·// 줄, 공백, 숫자 제거 후 대문자
  function clean(raw) {
    return String(raw || '').split(/\r?\n/)
      .filter(function (l) { var x = l.trim(); return x.charAt(0) !== '>' && !/^(ORIGIN|\/\/)/i.test(x); })
      .join('').replace(/[\s0-9]/g, '').toUpperCase();
  }
  function check(s) {
    if (!s) return '서열이 비어 있어요.';
    if (/U/.test(s)) return 'U가 들어 있어요. RNA 서열이라면 U를 T로 바꿔서 넣어 주세요.';
    var bad = s.replace(/[ACGTNRYSWKMBDHV]/g, '');
    if (bad) return 'DNA가 아닌 글자가 있어요: ' + Array.from(new Set(bad.split(''))).join(' ');
    return '';
  }
  function comp(s) { return s.split('').map(function (b) { return COMP[b] || 'N'; }).join(''); }
  function rev(s) { return s.split('').reverse().join(''); }
  function rc(s) { return rev(comp(s)); }
  function count(s, re) { return (s.match(re) || []).length; }
  function gcPercent(s) { return s.length ? count(s, /[GC]/g) / s.length * 100 : 0; }

  // ── Tm 계산 ─────────────────────────────────────────────
  function tmWallace(s) { return 2 * count(s, /[AT]/g) + 4 * count(s, /[GC]/g); }
  function tmBasic(s) { return 64.9 + 41 * (count(s, /[GC]/g) - 16.4) / s.length; }

  // SantaLucia (1998) unified nearest-neighbor 값. [ΔH kcal/mol, ΔS cal/(K·mol)]
  var NN = {
    AA: [-7.9, -22.2], TT: [-7.9, -22.2],
    AT: [-7.2, -20.4], TA: [-7.2, -21.3],
    CA: [-8.5, -22.7], TG: [-8.5, -22.7],
    GT: [-8.4, -22.4], AC: [-8.4, -22.4],
    CT: [-7.8, -21.0], AG: [-7.8, -21.0],
    GA: [-8.2, -22.2], TC: [-8.2, -22.2],
    CG: [-10.6, -27.2], GC: [-9.8, -24.4],
    GG: [-8.0, -19.9], CC: [-8.0, -19.9]
  };
  var NN_LABEL = {
    AA: 'AA/TT', TT: 'AA/TT', AT: 'AT/TA', TA: 'TA/AT', CA: 'CA/GT', TG: 'CA/GT', GT: 'GT/CA', AC: 'GT/CA',
    CT: 'CT/GA', AG: 'CT/GA', GA: 'GA/CT', TC: 'GA/CT', CG: 'CG/GC', GC: 'GC/CG', GG: 'GG/CC', CC: 'GG/CC'
  };

  // 기본 반응 조건: 일반적인 Taq PCR 버퍼(KCl 50 mM, MgCl2 1.5 mM, dNTP 각 0.2 mM, 프라이머 0.5 µM)
  var TAQ = { mono: 50, mg: 1.5, dntp: 0.8, primer: 0.5 };
  // 설계 기준 조건: 프라이머 설계 프로그램들이 Tm을 계산할 때 흔히 쓰는 값(프라이머 50 nM)
  var DESIGN = { mono: 50, mg: 1.5, dntp: 0.6, primer: 0.05 };

  function tmNN(s, cond) {
    cond = cond || TAQ;
    if (/[^ACGT]/.test(s) || s.length < 2) return null;
    var dH = 0, dS = 0, steps = [];
    for (var i = 0; i < s.length - 1; i++) {
      var p = NN[s.substr(i, 2)];
      dH += p[0]; dS += p[1];
      steps.push(s.substr(i, 2));
    }
    [s.charAt(0), s.charAt(s.length - 1)].forEach(function (b) {
      if (b === 'G' || b === 'C') { dH += 0.1; dS += -2.8; } else { dH += 2.3; dS += 4.1; }
    });
    var selfComp = s === rc(s);
    if (selfComp) dS += -1.4;
    // 2가 양이온 보정(von Ahsen 2001): Mg2+ 중 dNTP에 묶이지 않은 양을 1가 이온으로 환산
    var freeMg = Math.max(0, cond.mg - cond.dntp);
    var naEq = cond.mono + 120 * Math.sqrt(freeMg);           // mM
    var dSsalt = dS + 0.368 * (s.length - 1) * Math.log(naEq / 1000);
    var C = cond.primer * 1e-6;                                // M
    var R = 1.987;
    var tm = dH * 1000 / (dSsalt + R * Math.log(selfComp ? C : C / 4)) - 273.15;
    return { tm: tm, dH: dH, dS: dS, dSsalt: dSsalt, naEq: naEq, steps: steps, selfComp: selfComp };
  }

  // ── 2차 구조 · 반복 점검(간이) ─────────────────────────
  function pairs(a, b) { return COMP[a] === b; }

  // a(5'→3')와 b(5'→3')가 반대 방향으로 붙을 때 가장 긴 연속 짝과 3′ 끝이 걸린 짝
  function dimer(a, b) {
    var br = rev(b), best = { run: 0 }, best3 = { run: 0 };
    for (var o = -(br.length - 1); o < a.length; o++) {
      var run = 0, start = 0;
      for (var i = Math.max(0, o); i < Math.min(a.length, o + br.length); i++) {
        var j = i - o;
        if (pairs(a[i], br[j])) {
          if (run === 0) start = i;
          run++;
          var info = { run: run, aStart: start, aEnd: i, offset: o };
          if (run > best.run) best = info;
          var touches3 = i === a.length - 1 || (j - run + 1) === 0;
          if (touches3 && run > best3.run) best3 = info;
        } else run = 0;
      }
    }
    return { best: best, best3: best3 };
  }

  // 머리핀: 프라이머 한 가닥 안에서 접혀 붙는 줄기(stem). 고리(loop)는 3 nt 이상
  function hairpin(s) {
    var L = s.length, best = { stem: 0 };
    for (var i = 0; i < L; i++) {
      for (var j = L - 1; j > i; j--) {
        // 고리가 3 nt 미만이 되기 직전까지만 짝을 늘린다(끝까지 늘린 뒤 버리면 짧은 줄기를 놓침)
        var k = 0;
        while (j - i - 2 * k - 1 >= 3 && pairs(s[i + k], s[j - k])) k++;
        var loop = j - i - 2 * k + 1;
        if (k > best.stem) best = { stem: k, i: i, j: j, loop: loop, end3: j === L - 1 };
      }
    }
    return best;
  }

  function palindromes(s, min) {
    var out = [];
    for (var len = 12; len >= (min || 6); len -= 2) {
      for (var i = 0; i + len <= s.length; i++) {
        var w = s.substr(i, len);
        if (w === rc(w) && !out.some(function (p) { return i >= p.pos && i + len <= p.pos + p.seq.length; })) out.push({ pos: i, seq: w });
      }
    }
    return out;
  }

  // 자주 쓰는 제한효소(인식 서열, 자르는 위치: 인식 서열 첫 염기 앞을 0으로)
  var ENZYMES = [
    { name: 'EcoRI', site: 'GAATTC', cut: 1 }, { name: 'SacI', site: 'GAGCTC', cut: 5 },
    { name: 'KpnI', site: 'GGTACC', cut: 5 }, { name: 'SmaI', site: 'CCCGGG', cut: 3 },
    { name: 'BamHI', site: 'GGATCC', cut: 1 }, { name: 'XbaI', site: 'TCTAGA', cut: 1 },
    { name: 'SalI', site: 'GTCGAC', cut: 1 }, { name: 'PstI', site: 'CTGCAG', cut: 5 },
    { name: 'SphI', site: 'GCATGC', cut: 5 }, { name: 'HindIII', site: 'AAGCTT', cut: 1 },
    { name: 'NdeI', site: 'CATATG', cut: 2 }, { name: 'XhoI', site: 'CTCGAG', cut: 1 },
    { name: 'NotI', site: 'GCGGCCGC', cut: 2 }, { name: 'NcoI', site: 'CCATGG', cut: 1 },
    { name: 'BglII', site: 'AGATCT', cut: 1 }, { name: 'EcoRV', site: 'GATATC', cut: 3 },
    { name: 'NheI', site: 'GCTAGC', cut: 1 }, { name: 'SpeI', site: 'ACTAGT', cut: 1 }
  ];
  function findSites(s, circular) {
    var t = circular ? s + s.substr(0, 11) : s, hits = [];
    ENZYMES.forEach(function (e) {
      var i = t.indexOf(e.site);
      while (i >= 0 && i < s.length) { hits.push({ name: e.name, site: e.site, pos: i + 1, cutAfter: ((i + e.cut - 1) % s.length) + 1, enzyme: e }); i = t.indexOf(e.site, i + 1); }
    });
    return hits.sort(function (a, b) { return a.pos - b.pos; });
  }

  // 한 프라이머 점검. level: ok / warn / bad
  function qc(s) {
    var items = [], L = s.length, gc = gcPercent(s);
    function add(level, title, detail) { items.push({ level: level, title: title, detail: detail }); }

    add(L >= 18 && L <= 25 ? 'ok' : (L >= 15 && L <= 30 ? 'warn' : 'bad'), '길이 ' + L + ' nt', '권장 18–25 nt');
    add(gc >= 40 && gc <= 60 ? 'ok' : (gc >= 30 && gc <= 70 ? 'warn' : 'bad'), 'GC 비율 ' + gc.toFixed(0) + '%', '권장 40–60%');

    var last5 = s.substr(-5), gc5 = count(last5, /[GC]/g), lastGC = /[GC]$/.test(s);
    if (gc5 > 3) add('warn', '3′ 끝 마지막 5개 중 G/C ' + gc5 + '개', '3′ 끝이 너무 단단하면 비슷한 다른 자리에도 붙어 연장될 수 있어요. 1–3개가 적당합니다.');
    else if (gc5 === 0) add('warn', '3′ 끝 마지막 5개가 모두 A/T', '3′ 끝이 약하게 붙어 연장이 잘 시작되지 않을 수 있어요.');
    else add('ok', '3′ 끝 GC ' + gc5 + '개' + (lastGC ? ', G/C로 끝남' : ''), '마지막 5개 중 1–3개');

    var run = s.match(/(A{4,}|C{4,}|G{4,}|T{4,})/);
    if (run) add(run[0].length >= 5 ? 'bad' : 'warn', '같은 염기 반복: ' + run[0], '같은 염기가 4개 이상 이어지면 프라이머가 미끄러져 엉뚱하게 붙을 수 있어요.');
    var di = s.match(/(AC|AG|AT|CA|CG|CT|GA|GC|GT|TA|TC|TG)\1{3,}/);
    if (di) add('warn', '두 염기 반복: ' + di[0], '같은 두 염기가 4번 이상 반복되면 미끄러짐이 생길 수 있어요.');
    if (!run && !di) add('ok', '반복 서열 없음', '');

    var hp = hairpin(s);
    if (hp.stem >= 4) add(hp.stem >= 6 || hp.end3 ? 'bad' : 'warn', '머리핀(hairpin) 줄기 ' + hp.stem + ' bp' + (hp.end3 ? ', 3′ 끝 포함' : ''), '프라이머가 스스로 접히면 주형에 붙지 못해요.');
    else add('ok', '머리핀 없음', '');

    var d = dimer(s, s);
    if (d.best3.run >= 4) add('bad', '자기끼리 3′ 끝이 ' + d.best3.run + ' bp 붙음(self 3′-dimer)', '3′ 끝끼리 붙으면 중합효소가 서로를 늘려 프라이머 이합체가 생겨요.');
    else if (d.best.run >= 6) add('warn', '자기끼리 ' + d.best.run + ' bp 연속으로 붙음(self-dimer)', '같은 프라이머 두 개가 서로 붙을 수 있어요.');
    else add('ok', '자기끼리 붙는 구간 ' + d.best.run + ' bp 이하', '');

    var pal = palindromes(s, 6);
    pal.forEach(function (p) {
      var enz = ENZYMES.filter(function (e) { return p.seq.indexOf(e.site) >= 0; }).map(function (e) { return e.name; });
      add('warn', '회문 서열 ' + p.seq + ' (' + (p.pos + 1) + '번째부터)' + (enz.length ? ' = ' + enz.join(', ') + ' 자리' : ''),
        '앞뒤로 읽어도 역상보가 같은 서열이라 자기 자신과 그대로 짝을 이뤄요. 일반 PCR 프라이머에서는 피하는 게 좋고, 클로닝용 프라이머에서는 일부러 넣습니다.');
    });
    return { items: items, dimer: d, hairpin: hp };
  }

  function dimerPicture(a, b, info) {
    // a는 5'→3', b는 아래에 3'→5'로 표시
    if (!info || !info.run) return '';
    var br = rev(b), o = info.offset, pad = Math.max(0, -o);
    var top = '5′-' + ' '.repeat(pad) + a + '-3′';
    var bot = '3′-' + ' '.repeat(Math.max(0, o)) + br + '-5′';
    var mid = '   ';
    var len = Math.max(a.length + pad, br.length + Math.max(0, o));
    for (var x = 0; x < len; x++) {
      var i = x - pad, j = x - Math.max(0, o);
      mid += (i >= info.aStart && i <= info.aEnd && j >= 0 && j < br.length && pairs(a[i], br[j])) ? '|' : ' ';
    }
    return top + '\n' + mid + '\n' + bot;
  }

  // Taq 기준 PCR 조건
  function taqProgram(tmF, tmR, size) {
    var ta = Math.round(Math.min(tmF, tmR) - 5);
    ta = Math.max(50, Math.min(68, ta));
    var ext = size <= 500 ? 30 : Math.ceil(size / 1000 * 60 / 15) * 15;
    return { ta: ta, annealSec: 30, extSec: ext };
  }

  // ── mRNA–gDNA 맞대기(엑손 찾기) ─────────────────────────
  // gDNA의 16-mer 위치표를 만들고, mRNA를 왼쪽부터 훑으며 같은 순서로 나오는 조각을 엑손으로 잡는다
  function exonAlign1(M, G) {
    var K = 16, idx = {}, blocks = [], gPrev = -1, mPrev = -1, i = 0;
    for (var q = 0; q + K <= G.length; q++) { var w = G.substr(q, K); (idx[w] || (idx[w] = [])).push(q); }
    while (i + K <= M.length) {
      var cands = idx[M.substr(i, K)], g = -1;
      if (cands) for (var c = 0; c < cands.length; c++) if (cands[c] > gPrev) { g = cands[c]; break; }
      if (g < 0) { i++; continue; }
      var s = i, gs = g;
      while (s - 1 > mPrev && gs - 1 > gPrev && M[s - 1] === G[gs - 1]) { s--; gs--; }
      var len = K + (i - s);
      while (s + len < M.length && gs + len < G.length && M[s + len] === G[gs + len]) len++;
      blocks.push({ ms: s, me: s + len - 1, gs: gs, ge: gs + len - 1 });
      mPrev = s + len - 1; gPrev = gs + len - 1; i = s + len;
    }
    var ex = [];
    blocks.forEach(function (b) {
      var p = ex[ex.length - 1];
      if (p && (b.gs - b.ms) === (p.gs - p.ms) && b.ms - p.me - 1 <= 12) { p.me = b.me; p.ge = b.ge; }
      else ex.push({ ms: b.ms, me: b.me, gs: b.gs, ge: b.ge });
    });
    ex = ex.filter(function (b) { return b.me - b.ms + 1 >= 18; });
    // 경계 바로 옆 염기가 양쪽 엑손 어디에나 맞을 때는 경계를 밀 수 있다.
    // 그 범위 안에서 인트론이 GT로 시작해 AG로 끝나는 자리(스플라이싱 규칙)로 경계를 옮긴다. NCBI 레코드도 이 규칙을 따른다.
    for (var k = 0; k + 1 < ex.length; k++) {
      var e1 = ex[k], e2 = ex[k + 1], o1 = e1.gs - e1.ms, o2 = e2.gs - e2.ms;
      if (e2.ms !== e1.me + 1) continue;
      var b0 = e1.me, best = null;
      for (var d = -12; d <= 12; d++) {
        var b = b0 + d, ok = b > e1.ms && b < e2.me;
        for (var x = Math.min(b, b0) + 1; ok && x <= Math.max(b, b0); x++) ok = d < 0 ? M[x] === G[x + o2] : M[x] === G[x + o1];
        if (!ok) continue;
        var donor = G.substr(b + o1 + 1, 2), acc = G.substr(b + 1 + o2 - 2, 2);
        var sc = (donor === 'GT' ? 2 : donor === 'GC' ? 1 : 0) + (acc === 'AG' ? 2 : 0);
        if (!best || sc > best.sc || (sc === best.sc && Math.abs(d) < Math.abs(best.d))) best = { sc: sc, d: d };
      }
      if (best && best.d && best.sc >= 3) {
        e1.me += best.d; e1.ge += best.d; e2.ms += best.d; e2.gs += best.d;
      }
    }
    var covered = ex.reduce(function (a, b) { return a + (b.me - b.ms + 1); }, 0);
    return { exons: ex, cover: M.length ? covered / M.length : 0, G: G, flipped: false };
  }
  function exonAlign(M, G) {
    var A = exonAlign1(M, G);
    if (A.cover < 0.5) { var B = exonAlign1(M, rc(G)); if (B.cover > A.cover) { B.flipped = true; A = B; } }
    return A;
  }
  function exonOf(ex, pos) { for (var i = 0; i < ex.length; i++) if (pos >= ex[i].ms && pos <= ex[i].me) return i; return -1; }

  var CODON = {};
  (function () {
    var b = 'TCAG', aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', n = 0;
    for (var x = 0; x < 4; x++) for (var y = 0; y < 4; y++) for (var z = 0; z < 4; z++) CODON[b[x] + b[y] + b[z]] = aa[n++];
  })();
  var AA3 = { A: 'Ala', R: 'Arg', N: 'Asn', D: 'Asp', C: 'Cys', Q: 'Gln', E: 'Glu', G: 'Gly', H: 'His', I: 'Ile', L: 'Leu', K: 'Lys', M: 'Met', F: 'Phe', P: 'Pro', S: 'Ser', T: 'Thr', W: 'Trp', Y: 'Tyr', V: 'Val', '*': '종결' };
  function translate(s) { var p = ''; for (var i = 0; i + 3 <= s.length; i += 3) p += CODON[s.substr(i, 3)] || 'X'; return p; }

  root.DNA = {
    clean: clean, check: check, comp: comp, rev: rev, rc: rc, gcPercent: gcPercent,
    tmWallace: tmWallace, tmBasic: tmBasic, tmNN: tmNN, NN: NN, NN_LABEL: NN_LABEL, TAQ: TAQ, DESIGN: DESIGN,
    dimer: dimer, hairpin: hairpin, palindromes: palindromes, qc: qc, dimerPicture: dimerPicture,
    ENZYMES: ENZYMES, findSites: findSites, taqProgram: taqProgram, CODON: CODON, AA3: AA3, translate: translate,
    exonAlign: exonAlign, exonOf: exonOf
  };
})(window);
