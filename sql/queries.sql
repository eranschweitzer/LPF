SELECT v.casename, v.id, v.criteria, v.prop, v.value vval, a.prop, a.value aval , pf.prop, pf.value pfval, pt.prop, pt.value ptval, qf.prop, qf.value qfval, qt.prop, qt.value qtval
FROM linpf.bestnogw1 v, linpf.bestnogw1 a, linpf.bestnogw1 pf, linpf.bestnogw1 pt, linpf.bestnogw1 qf, linpf.bestnogw1 qt
WHERE v.casename=a.casename and v.casename=pf.casename and v.casename=pt.casename and v.casename=qf.casename and v.casename=qt.casename 
and v.id=a.id and v.id=pf.id and v.id=pt.id and v.id=qf.id and v.id=qt.id 
and v.criteria='max' and v.criteria=a.criteria and v.criteria=pf.criteria and v.criteria=pt.criteria and v.criteria=qf.criteria and v.criteria=qt.criteria 
and v.prop='volt' and a.prop='ang' and pf.prop like 'pf%' and pt.prop like 'pt%' and qf.prop like 'qf%' and qt.prop like 'qt%'
ORDER BY v.casename
;

SELECT v.casename, v.id, v.criteria, v.prop, v.value vval, a.prop, a.value aval
FROM linpf.bestnogw1 v, linpf.bestnogw1 a
WHERE v.casename=a.casename 
and v.id=a.id 
and v.criteria='max' and v.criteria=a.criteria 
and v.prop='volt' and a.prop='ang' 
ORDER BY v.casename
;

SELECT pf.casename, pf.id, pf.criteria, pf.prop, pf.value pfval, qf.prop, qf.value qfval 
FROM linpf.bestnogw1 pf, linpf.bestnogw1 qf
WHERE pf.casename=qf.casename
and pf.id=qf.id
and pf.criteria='max' and pf.criteria=qf.criteria
and pf.prop like 'pf%' and qf.prop like 'qf%'
ORDER BY pf.casename
;

SELECT pf.casename, pf.id, pf.criteria, pf.prop, pf.value pfval, qf.prop, qf.value qfval 
FROM linpf.bestnogw0 pf, linpf.bestnogw0 qf
WHERE pf.casename=qf.casename
and pf.id=qf.id
and pf.criteria='max' and pf.criteria=qf.criteria
and pf.prop like 'pf%' and qf.prop like 'qf%'
ORDER BY pf.casename
;

SELECT DISTINCT(casename), sum(dc::INTEGER) OVER(partition by casename) dc, sum(linpf::INTEGER) OVER(partition by casename) linpf FROM linpf.actestnogw ORDER BY dc ;

WITH tmp AS (
SELECT v.casename, v.id, v.criteria, v.prop, v.value, a.dc, a.linpf, ROW_NUMBER() OVER (PARTITION BY v.casename ORDER BY v.id) as rk 
FROM linpf.bestnogw1 v, linpf.actestnogw a
WHERE v.casename=a.casename and v.id=a.id and a.linpf is TRUE and v.criteria='max' and v.prop like 'pf%'
ORDER BY v.casename)
SELECT * FROM tmp WHERE rk=1;

WITH cor as (
		SELECT id, count(*) as cnt FROM linpf.bestnogw1 where prop='volt' and criteria='cor' GROUP BY id), 
	avg as (
		SELECT id, count(*) as cnt FROM linpf.bestnogw1 where prop='volt' and criteria='avg' GROUP BY id), 
	max as (
		SELECT id, count(*) as cnt FROM linpf.bestnogw1 where prop='volt' and criteria='max' GROUP BY id) 
	SELECT cor.id, cor.cnt as corr_cnt, avg.cnt as avg_cnt, max.cnt as max_cnt FROM cor, avg, max WHERE cor.id = avg.id and avg.id = max.id;

WITH cor as (
		SELECT id, count(*) as cnt FROM linpf.bestnogw1 where prop='ang' and criteria='cor' GROUP BY id), 
	avg as (
		SELECT id, count(*) as cnt FROM linpf.bestnogw1 where prop='ang' and criteria='avg' GROUP BY id), 
	max as (
		SELECT id, count(*) as cnt FROM linpf.bestnogw1 where prop='ang' and criteria='max' GROUP BY id) 
	SELECT cor.id, cor.cnt as corr_cnt, avg.cnt as avg_cnt, max.cnt as max_cnt FROM cor, avg, max WHERE cor.id = avg.id and avg.id = max.id

WITH vlt as (
		SELECT id, count(*) as cnt FROM linpf.bestnogw1 where prop='volt' and criteria='max' GROUP BY id), 
	ang as (
		SELECT id, count(*) as cnt FROM linpf.bestnogw1 where prop='ang' and criteria='max' GROUP BY id) 
	SELECT vlt.id, vlt.cnt as vlt_cnt, ang.cnt as ang_cnt FROM vlt, ang WHERE vlt.id = ang.id;

-- ====================================================================
-- From new formulation (tables are linpf.vt, linpf.flows, linpf.cases)
-- ====================================================================

WITH tmp AS (
	SELECT caseid, criteria, min(value) as v FROM linpf.vt WHERE prop='vol' GROUP BY caseid, criteria ORDER BY caseid)
SELECT b.name, a.criteria, a.id FROM linpf.vt a, tmp, linpf.cases b WHERE b.id=a.caseid and a.value=tmp.v and a.criteria=tmp.criteria ORDER BY b.name;
