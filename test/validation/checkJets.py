"""Checks the NANO output of validateJets_cfg.py. Usage: python3 -I checkJets.py jets_nano.root"""
import sys
import numpy as np
import uproot
import awkward as ak

t = uproot.open(sys.argv[1])['Events']
a = t.arrays(library='ak')
ok = True


def report(name, cond, detail):
    global ok
    ok &= bool(cond)
    print(('PASS ' if cond else 'FAIL ') + name + ': ' + detail)


def sums(prefix):
    c = prefix + 'Constituents_'
    pt, phi, eta, ij = a[c + 'pt'], a[c + 'phi'], a[c + 'eta'], a[c + 'iJet']
    out = []
    for ev in range(len(pt)):
        n = len(a[prefix + '_pt'][ev])
        px = np.bincount(np.asarray(ij[ev]), weights=np.asarray(pt[ev] * np.cos(phi[ev])), minlength=n)
        py = np.bincount(np.asarray(ij[ev]), weights=np.asarray(pt[ev] * np.sin(phi[ev])), minlength=n)
        out.append(np.hypot(px, py))
    return ak.Array(out)


def flat(x):
    return np.asarray(ak.flatten(x))


# (a) reco PAT jets: PUPPI-weighted constituent sum reproduces the raw jet pt
r = flat(sums('RecoAK8')) / flat(a.RecoAK8_rawPt)
report('reco AK8 constituent sum / rawPt', np.all(np.abs(r - 1) < 0.01) and len(r) >= 10,
       'n=%d median %.5f max|dev| %.4f' % (len(r), np.median(r), np.max(np.abs(r - 1))))
# (b) gen jets: constituent sum reproduces the gen-jet pt
g = flat(sums('GenAK8')) / flat(a.GenAK8_pt)
# MiniAOD gen-jet momenta differ from the sum of their packed daughters by up to ~3e-4 (packing precision).
report('gen AK8 constituent sum / pt', np.all(np.abs(g - 1) < 1e-3), 'n=%d max|dev| %.2e' % (len(g), np.max(np.abs(g - 1))))


def matchRatio(ref, other, refPt, minPt=200):
    ratios = []
    for ev in range(len(a)):
        for pt, eta, phi in zip(a[refPt][ev], a[ref + '_eta'][ev], a[ref + '_phi'][ev]):
            if pt < minPt:
                continue
            dr = np.hypot(np.asarray(a[other + '_eta'][ev]) - eta,
                          np.remainder(np.asarray(a[other + '_phi'][ev]) - phi + np.pi, 2 * np.pi) - np.pi)
            if len(dr) and dr.min() < 0.05:
                ratios.append(a[other + '_pt'][ev][int(dr.argmin())] / pt)
            else:
                ratios.append(np.nan)
    return np.array(ratios)


# (c) reclustered gen AK8 (no neutrinos) reproduces slimmedGenJetsAK8
m = matchRatio('GenAK8', 'GenAK8Rc', 'GenAK8_pt')
report('gen recluster vs slimmedGenJetsAK8', np.nanmax(np.abs(m - 1)) < 1e-3 and np.mean(np.isnan(m)) < 0.01,
       'n=%d unmatched %d max|dev| %.2e' % (len(m), np.isnan(m).sum(), np.nanmax(np.abs(m - 1))))
# (d) reclustered PUPPI AK8 vs slimmedJetsAK8 raw pt
m = matchRatio('RecoAK8', 'RecoAK8Rc', 'RecoAK8_rawPt')
report('reco recluster vs slimmedJetsAK8 raw', np.mean(np.isnan(m)) < 0.05 and np.nanmedian(np.abs(m - 1)) < 0.01,
       'n=%d unmatched %d median|dev| %.2e p95|dev| %.2e' % (len(m), np.isnan(m).sum(), np.nanmedian(np.abs(m - 1)),
                                                             np.nanpercentile(np.abs(m - 1), 95)))

# variations: same source jets before the jet selection may differ, so compare by srcIdx
def bySrc(prefix, col):
    return [dict(zip(np.asarray(s).tolist(), np.asarray(v).tolist())) for s, v in zip(a[prefix + '_srcIdx'], a[prefix + '_' + col])]


nom, up = bySrc('RecoAK8', 'pt'), bySrc('RecoAK8JESUp', 'pt')
rj = np.array([u[k] / n[k] for n, u in zip(nom, up) for k in n if k in u])
report('JESUp / nominal', np.all(rj > 1) and np.all(rj < 1.1), 'n=%d range %.4f-%.4f' % (len(rj), rj.min(), rj.max()))
jup = bySrc('RecoAK8JERUp', 'pt')
rr = np.array([u[k] / n[k] for n, u in zip(nom, jup) for k in n if k in u])
report('JERUp differs from nominal', np.mean(np.abs(rr - 1) > 1e-6) > 0.5, 'n=%d frac changed %.2f' % (len(rr), np.mean(np.abs(rr - 1) > 1e-6)))
# raw vs nominal (JER nominal applied): changed but small
nr = np.array([n[k] / raw[k] for n, raw in zip(nom, bySrc('RecoAK8', 'rawPt')) for k in n])
print('INFO  nominal pt / rawPt (JEC x JER): median %.3f' % np.median(nr))

cq = flat(a.RecoAK8Constituents_charge); cpdg = np.abs(flat(a.RecoAK8Constituents_pdgId))
nCharged = np.sum((cq != 0) & (cpdg != 11) & (cpdg != 13))
dq = flat(a.RecoAK8DropChargedConstituents_charge); dpdg = np.abs(flat(a.RecoAK8DropChargedConstituents_pdgId))
dCharged = np.sum((dq != 0) & (dpdg != 11) & (dpdg != 13))
nNeutral = len(cq) - nCharged
dNeutral = len(dq) - dCharged
# The JER-nominal jets feed both chains, but the jet selection (pt > 0) keeps all jets in both, so counts compare.
report('DropCharged keeps ~half the charged hadrons and all others', abs(dCharged / nCharged - 0.5) < 0.02 and dNeutral == nNeutral,
       'charged kept %.3f, others %d vs %d' % (dCharged / nCharged, dNeutral, nNeutral))
ph = flat(a.RecoAK8Constituents_pt)[cpdg == 22]
sph = flat(a.RecoAK8ScalePhotonConstituents_pt)[np.abs(flat(a.RecoAK8ScalePhotonConstituents_pdgId)) == 22]
report('ScalePhoton scales photon pt by 1.1', len(ph) == len(sph) and np.allclose(sph / ph, 1.1, rtol=1e-5),
       'n=%d' % len(ph))
ids = flat(a.RecoAK8_jetId)
tid = flat(a.RecoAK8TightId_jetId)
report('TightId keeps only jetId bit 1', np.all((tid & 2) == 2) and len(tid) == np.sum((ids & 2) == 2),
       'nominal %d, tightLepVeto %d, tight %d' % (len(ids), len(tid), np.sum((ids & 1) == 1)))
sys.exit(0 if ok else 1)
