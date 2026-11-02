/* 11주차 도메인 교환 클로닝 · 공통 계산
 *
 * 모든 계산은 브라우저 안에서만 이루어지고, 서열은 어디에도 전송되지 않는다.
 * dna.js(10주차와 같은 파일)와 receptors.js를 먼저 읽어야 한다.
 */
(function (root) {
  var D = root.DNA, R = root.RECEPTORS, T = root.TASK;

  /* ── 제한효소 ──────────────────────────────────────────
   * cut: 인식 서열 첫 염기 앞을 0으로 둔, 위쪽 가닥을 자르는 위치
   * ends: 생기는 말단 모양. temp: 반응 온도 (공급사 자료 기준)
   */
  var ENZ = [
    { name: 'EcoRI',   site: 'GAATTC',   cut: 1, ends: "5′ 점착", temp: 37 },
    { name: 'SacI',    site: 'GAGCTC',   cut: 5, ends: "5′ 점착", temp: 37 },
    { name: 'KpnI',    site: 'GGTACC',   cut: 5, ends: "5′ 점착", temp: 37 },
    { name: 'XmaI',    site: 'CCCGGG',   cut: 1, ends: "5′ 점착", temp: 37 },
    { name: 'SmaI',    site: 'CCCGGG',   cut: 3, ends: '평활',    temp: 25 },
    { name: 'BamHI',   site: 'GGATCC',   cut: 1, ends: "5′ 점착", temp: 37 },
    { name: 'XbaI',    site: 'TCTAGA',   cut: 1, ends: "5′ 점착", temp: 37 },
    { name: 'SalI',    site: 'GTCGAC',   cut: 1, ends: "5′ 점착", temp: 37 },
    { name: 'PstI',    site: 'CTGCAG',   cut: 5, ends: "3′ 점착", temp: 37 },
    { name: 'SphI',    site: 'GCATGC',   cut: 5, ends: "3′ 점착", temp: 37 },
    { name: 'HindIII', site: 'AAGCTT',   cut: 1, ends: "5′ 점착", temp: 37 }
  ];
  function enz(name) { for (var i = 0; i < ENZ.length; i++) if (ENZ[i].name === name) return ENZ[i]; return null; }

  /* 두 효소를 한 통에 같이 넣을 수 있는지 */
  function canCombine(a, b) {
    if (!a || !b) return { ok: false, why: '효소를 두 개 고르세요.' };
    if (a.name === b.name) return { ok: false, why: '같은 효소를 양쪽에 쓰면 조각이 뒤집혀 들어갈 수 있습니다.' };
    if (a.temp !== b.temp) {
      return { ok: false, kind: 'temp',
        why: a.name + '은(는) ' + a.temp + ' °C, ' + b.name + '은(는) ' + b.temp + ' °C에서 일합니다. 한 반응에서 두 온도를 맞출 수 없습니다.' };
    }
    return { ok: true, why: '둘 다 ' + a.temp + ' °C에서 일하므로 한 반응에서 함께 처리할 수 있습니다.' };
  }

  /* 자른 자리 그림: GAATTC → G^AATTC / CTTAA^G */
  function cutPicture(e) {
    var top = e.site, bot = D.comp(top), cutTop = e.cut, cutBot = top.length - e.cut;
    return '5′- ' + top.slice(0, cutTop) + '▼' + top.slice(cutTop) + ' -3′\n' +
           '3′- ' + bot.slice(0, cutBot) + '▲' + bot.slice(cutBot) + ' -5′';
  }

  /* ── 수용체 서열에서 도메인 꺼내기 ───────────────────── */
  function cds(sym) { return R[sym].cds; }
  function aaRange(sym, from, to) { return cds(sym).substring((from - 1) * 3, to * 3); }
  function domainSeq(sym, key) {
    var d = R[sym].domains.filter(function (x) { return x.key === key; })[0];
    return d ? aaRange(sym, d.aaFrom, d.aaTo) : '';
  }
  /* 오늘 쓸 두 조각. B는 원래 CDS에 있는 종결코돈까지 가져온다 */
  function fragSeq(which) {
    var t = T[which], s = aaRange(t.sym, t.aaFrom, t.aaTo);
    if (t.addStop) {
      var next = cds(t.sym).substr(t.aaTo * 3, 3);
      s += (next.length === 3 && D.CODON[next] === '*') ? next : 'TGA';
    }
    return s;
  }
  function chimera() { return fragSeq('A') + fragSeq('B'); }

  /* Tm을 숫자로. dna.js의 tmNN은 계산 과정까지 담은 객체를 돌려준다 */
  function tm(seq, cond) {
    var s = String(seq || '').replace(/[^ACGT]/g, '');
    if (s.length < 8) return null;
    var r = D.tmNN(s, cond || D.TAQ);
    return r ? r.tm : null;
  }

  /* 번역틀 점검 */
  function frameCheck(seq) {
    var out = { len: seq.length, mod: seq.length % 3, stops: [], prot: D.translate(seq) };
    var p = out.prot;
    for (var i = 0; i < p.length; i++) if (p[i] === '*' && i < p.length - 1) out.stops.push(i + 1);
    out.ok = out.mod === 0 && out.stops.length === 0;
    return out;
  }

  /* ── pUC19 ──────────────────────────────────────────── */
  var MCS = { from: 396, to: 454 };
  function vecSites() {
    var S = root.PUC19.seq, L = S.length, out = [];
    ENZ.forEach(function (e) {
      var hits = [], i = (S + S.substr(0, 10)).indexOf(e.site);
      while (i >= 0 && i < L) { hits.push(i + 1); i = (S + S.substr(0, 10)).indexOf(e.site, i + 1); }
      out.push({ enzyme: e, name: e.name, hits: hits, n: hits.length,
        pos: hits.length ? hits[0] : null, inMCS: hits.length === 1 && hits[0] >= MCS.from && hits[0] <= MCS.to,
        cutAfter: hits.length ? hits[0] - 1 + e.cut : null });
    });
    return out;
  }
  /* MCS 안에서 한 번만 자르는 효소, 위치 순서대로 */
  function mcsEnzymes() {
    return vecSites().filter(function (v) { return v.inMCS; }).sort(function (a, b) { return a.pos - b.pos; });
  }

  /* ── 프라이머가 붙는 자리 찾기 ───────────────────────
   * 3′ 끝에서부터 SEED nt가 정확히 맞는 곳만 자리로 본다.
   * 꼬리(5′쪽 안 맞는 부분)는 그대로 둔다 — 꼬리 달린 프라이머가 바로 이렇게 붙는다.
   */
  var SEED = 12;
  function bindSites(primer, strand) {
    var p = primer.toUpperCase(), out = [];
    if (p.length < SEED) return out;
    var seed = p.slice(-SEED), i = strand.indexOf(seed);
    while (i >= 0) {
      var endIdx = i + SEED - 1;                 // 프라이머 3′ 끝이 닿는 위치(0부터)
      var m = SEED;                              // 3′에서 이어지는 일치 길이
      while (m < p.length) {
        var tj = endIdx - m, pj = p.length - 1 - m;
        if (tj < 0 || strand[tj] !== p[pj]) break;
        m++;
      }
      var mism = 0;
      for (var k = m; k < p.length; k++) {
        var tk = endIdx - k;
        if (tk >= 0 && strand[tk] !== p[p.length - 1 - k]) mism++;
      }
      out.push({ end: endIdx + 1, start: endIdx - p.length + 2, match: m, mism: mism,
        matched: p.slice(p.length - m), tail: p.slice(0, p.length - m) });
      i = strand.indexOf(seed, i + 1);
    }
    return out;
  }
  /* 붙는 세기: 맞는 부분의 Tm과 결합 온도를 비교 */
  function bindEff(site, ta) {
    var seq = site.matched.replace(/[^ACGT]/g, '');
    if (seq.length < SEED) return 0;
    var t = tm(seq);
    if (t == null) return 0;
    var e = 1 / (1 + Math.exp((ta - t + 2) / 1.6));
    return Math.max(0, Math.min(1, e));
  }

  /* ── PCR 한 판 ─────────────────────────────────────────
   * tmpl: [{name, seq, copies}] 통 안에 있는 주형들
   * prog: {ta, extSec, cycles}
   */
  var AVO = 6.022e23, RATE = 1000;          // Taq 연장 속도 bp/분
  /* 통 안에 처음 들어 있는 분자 수. 사이클을 줄이면 띠가 흐려지도록 잡은 값 */
  var COPIES = { cdna: 2e5, fused: 1e8, leftover: 1e9 };
  function ngOf(copies, bp) { return copies * bp * 650 / AVO * 1e9; }

  function pcr(tmpls, fwd, rev, prog) {
    var notes = [], prods = [];
    fwd = D.clean(fwd); rev = D.clean(rev);
    if (!fwd || !rev) return { bands: [], notes: ['프라이머를 두 개 다 넣어 주세요.'] };

    tmpls.forEach(function (t) {
      var plus = t.seq, minus = D.rc(t.seq);
      var fOn = bindSites(fwd, plus), rOn = bindSites(rev, minus);
      var fOnMinus = bindSites(fwd, minus), rOnPlus = bindSites(rev, plus);
      if (!fOn.length && fOnMinus.length) notes.push(t.name + ': Forward 프라이머가 아래쪽 가닥에 붙습니다. 윗가닥 서열을 그대로 쓰는지 확인하세요.');
      if (!rOn.length && rOnPlus.length) notes.push(t.name + ': Reverse 프라이머가 윗가닥에 붙습니다. 역상보로 바꿔야 합니다.');

      fOn.forEach(function (f) {
        rOn.forEach(function (r) {
          // minus 가닥 좌표를 plus 좌표로 바꾼다
          var rEndPlus = t.seq.length - r.end + 1;      // reverse 3′ 끝이 닿는 plus 위치
          var rStartPlus = t.seq.length - r.start + 1;  // reverse 5′ 끝 쪽
          if (rEndPlus < f.end) return;                 // 서로 마주 보지 않음
          var inner = t.seq.substring(f.end, rEndPlus - 1);
          var seq = fwd + inner + D.rc(rev);
          var bp = seq.length;
          if (bp > 12000) return;
          var eF = bindEff(f, prog.ta), eR = bindEff(r, prog.ta);
          var cap = RATE * prog.extSec / 60;
          var sizeF = 1 / (1 + Math.pow(bp / Math.max(60, cap), 6));
          var E = eF * eR * sizeF * 0.92;
          if (E < 0.02) return;
          var copies = t.copies * Math.pow(1 + E, prog.cycles);
          prods.push({ bp: bp, seq: seq, copies: copies, tmpl: t.name,
            eF: eF, eR: eR, sizeF: sizeF, fSite: f, rSite: r, rEndPlus: rEndPlus, fEnd: f.end });
        });
      });
    });

    // 프라이머 이합체: 3′ 끝이 걸려 있으면 서로를 늘려 작은 띠가 생긴다
    var run3 = D.dimer(fwd, rev).best3.run || 0;
    if (run3 >= 4) {
      var dimerBp = Math.max(30, fwd.length + rev.length - run3);
      prods.push({ bp: dimerBp, seq: '', copies: 1e9 * Math.pow(2.2, Math.min(8, run3 - 3)), tmpl: '프라이머 이합체', dimer: true });
      notes.push('두 프라이머의 3′ 끝이 ' + run3 + ' 염기 서로 결합합니다(이합체). 젤 아래쪽에 작은 띠가 생깁니다.');
    }

    // 총량 한계: dNTP가 떨어지면 더 안 늘어난다 (50 µL 반응에서 약 450 ng)
    var totalBp = prods.reduce(function (s, p) { return s + p.copies * p.bp; }, 0);
    var MAXBP = 4.2e14;
    var scale = totalBp > MAXBP ? MAXBP / totalBp : 1;
    prods.forEach(function (p) { p.ng = ngOf(p.copies * scale, p.bp); });

    // 같은 크기끼리 합치기
    prods.sort(function (a, b) { return b.ng - a.ng; });
    var bands = [];
    prods.forEach(function (p) {
      var hit = bands.filter(function (b) { return Math.abs(b.bp - p.bp) <= Math.max(6, b.bp * 0.01); })[0];
      if (hit) { hit.ng += p.ng; if (p.ng > hit.ng / 2 && !hit.dimer) hit.parts.push(p); }
      else bands.push({ bp: p.bp, ng: p.ng, seq: p.seq, dimer: !!p.dimer, tmpl: p.tmpl, parts: [p], eF: p.eF, eR: p.eR, sizeF: p.sizeF });
    });
    bands = bands.filter(function (b) { return b.ng >= 0.4; }).sort(function (a, b) { return b.bp - a.bp; });

    if (!bands.length) notes.push('띠가 하나도 보이지 않습니다. 프라이머가 주형에 붙지 못했거나, 결합 온도가 너무 높습니다.');
    prods.forEach(function (p) {
      if (p.sizeF < 0.5 && !p.dimer && p.ng < 8) notes.push(p.bp.toLocaleString() + ' bp를 만들기에는 연장 시간이 짧아요. ' + Math.ceil(p.bp / RATE * 60 / 15) * 15 + '초 이상 주세요.');
    });
    return { bands: bands, notes: Array.from(new Set(notes)) };
  }

  /* ── 젤 ────────────────────────────────────────────────
   * 1 kb 래더 기준. 로그 눈금으로 이동 거리를 계산한다.
   */
  var LADDER = [10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1000, 500, 250];
  var GEL = { top: 24, bottom: 400, hi: 12000, lo: 150 };
  function yOf(bp) {
    var b = Math.max(GEL.lo, Math.min(GEL.hi, bp));
    var f = (Math.log10(GEL.hi) - Math.log10(b)) / (Math.log10(GEL.hi) - Math.log10(GEL.lo));
    return GEL.top + f * (GEL.bottom - GEL.top);
  }
  function bpOf(y) {
    var f = (y - GEL.top) / (GEL.bottom - GEL.top);
    return Math.round(Math.pow(10, Math.log10(GEL.hi) - f * (Math.log10(GEL.hi) - Math.log10(GEL.lo))));
  }
  var VISIBLE_NG = 1.5;   // 젤에서 눈에 보이기 시작하는 양

  /* ── 자연 냉각 결합 + 한 번의 연장 (sewing 1단계) ────
   * 두 조각이 겹치는 부분으로 서로 붙은 뒤, 서로를 주형 삼아 끝까지 채워진다.
   */
  function overlapLen(a, b) {
    var max = Math.min(a.length, b.length, 60);
    for (var n = max; n >= 8; n--) if (a.slice(-n) === b.slice(0, n)) return n;
    return 0;
  }
  function anneal(a, b, opt) {
    opt = opt || {};
    var ov = overlapLen(a, b), out = { overlap: ov, notes: [] };
    if (!ov) {
      out.ok = false;
      out.notes.push('두 조각에 겹치는 서열이 없습니다. 이음새 쪽 프라이머의 5′에 상대 조각의 서열을 꼬리로 붙여야 합니다.');
      return out;
    }
    var seq = a + b.slice(ov);
    out.ok = true; out.seq = seq; out.bp = seq.length;
    out.tm = tm(a.slice(-ov));
    if (ov < 15) out.notes.push('겹치는 부분이 ' + ov + ' nt뿐이에요. 15 nt 아래면 서로 붙어 있는 시간이 짧아 잘 이어지지 않습니다.');
    if (opt.holdSec != null && opt.holdSec < 60) out.notes.push('연장 시간이 ' + opt.holdSec + '초예요. ' + Math.ceil(seq.length / RATE * 60 / 30) * 30 + '초 이상 주면 끝까지 채워집니다.');
    if (opt.coolMin != null && opt.coolMin < 3) out.notes.push('식히는 시간이 ' + opt.coolMin + '분이면 너무 빨라요. 천천히 식혀야 겹친 부분끼리 제대로 만납니다.');
    return out;
  }

  /* ── 저장 (단계 사이로 결과 넘기기) ──────────────────── */
  var KEY = 'dsl11';
  function save(patch) {
    var cur = load();
    Object.keys(patch).forEach(function (k) { cur[k] = patch[k]; });
    try { localStorage.setItem(KEY, JSON.stringify(cur)); } catch (e) {}
    return cur;
  }
  function load() {
    try { return JSON.parse(localStorage.getItem(KEY) || '{}') || {}; } catch (e) { return {}; }
  }
  function clearAll() { try { localStorage.removeItem(KEY); } catch (e) {} }

  /* 상자 아래 전체를 끌어 높이를 바꾸는 손잡이 */
  function attachGrip(box, minH) {
    if (!box || box.__grip) return;
    box.__grip = true;
    var grip = document.createElement('div');
    grip.className = 'grip';
    grip.setAttribute('title', '끌어서 높이 바꾸기');
    box.parentNode.insertBefore(grip, box.nextSibling);
    var startY = 0, startH = 0, on = false;
    grip.addEventListener('pointerdown', function (ev) {
      on = true; startY = ev.clientY; startH = box.getBoundingClientRect().height;
      grip.classList.add('dragging');
      try { grip.setPointerCapture(ev.pointerId); } catch (e) {}
      ev.preventDefault();
    });
    grip.addEventListener('pointermove', function (ev) {
      if (!on) return;
      box.style.height = Math.max(minH || 120, startH + (ev.clientY - startY)) + 'px';
      ev.preventDefault();
    });
    function stop(ev) {
      if (!on) return;
      on = false; grip.classList.remove('dragging');
      try { grip.releasePointerCapture(ev.pointerId); } catch (e) {}
    }
    grip.addEventListener('pointerup', stop);
    grip.addEventListener('pointercancel', stop);
  }

  /* 숫자 보기 좋게 */
  function bpText(bp) { return bp >= 1000 ? (bp / 1000).toFixed(bp >= 10000 ? 1 : 2) + ' kb' : bp + ' bp'; }
  function ngText(ng) { return ng >= 100 ? Math.round(ng) + ' ng' : ng >= 10 ? ng.toFixed(1) + ' ng' : ng.toFixed(2) + ' ng'; }

  root.CLONE = {
    ENZ: ENZ, enz: enz, canCombine: canCombine, cutPicture: cutPicture,
    cds: cds, aaRange: aaRange, domainSeq: domainSeq, fragSeq: fragSeq, chimera: chimera, frameCheck: frameCheck,
    MCS: MCS, vecSites: vecSites, mcsEnzymes: mcsEnzymes,
    tm: tm, SEED: SEED, bindSites: bindSites, bindEff: bindEff, pcr: pcr, ngOf: ngOf,
    LADDER: LADDER, GEL: GEL, yOf: yOf, bpOf: bpOf, VISIBLE_NG: VISIBLE_NG, RATE: RATE, COPIES: COPIES,
    overlapLen: overlapLen, anneal: anneal,
    save: save, load: load, clearAll: clearAll, bpText: bpText, ngText: ngText, attachGrip: attachGrip
  };
})(window);
