CREATE TABLE linpf.best AS
WITH tmp as (
(SELECT casename, 'pf' as prop, criteria, min(value) as val FROM linpf.res WHERE prop like 'pf%' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'pf' as prop, criteria, max(value) as val FROM linpf.res WHERE prop like 'pf%' and criteria in ('cor') GROUP BY casename, criteria )
UNION
(SELECT casename, 'pt' as prop, criteria, min(value) as val FROM linpf.res WHERE prop like 'pt%' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'pt' as prop, criteria, max(value) as val FROM linpf.res WHERE prop like 'pt%' and criteria in ('cor') GROUP BY casename, criteria )
UNION
(SELECT casename, 'qf' as prop, criteria, min(value) as val FROM linpf.res WHERE prop like 'qf%' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'qf' as prop, criteria, max(value) as val FROM linpf.res WHERE prop like 'qf%' and criteria in ('cor') GROUP BY casename, criteria )
UNION
(SELECT casename, 'qt' as prop, criteria, min(value) as val FROM linpf.res WHERE prop like 'qt%' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'qt' as prop, criteria, max(value) as val FROM linpf.res WHERE prop like 'qt%' and criteria in ('cor') GROUP BY casename, criteria )
UNION
(SELECT casename, 'volt' as prop, criteria, min(value) as val FROM linpf.res WHERE prop ='volt' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'volt' as prop, criteria, max(value) as val FROM linpf.res WHERE prop ='volt' and criteria in ('cor') GROUP BY casename, criteria )
UNION
(SELECT casename, 'ang' as prop, criteria, min(value) as val FROM linpf.res WHERE prop ='ang' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'ang' as prop, criteria, max(value) as val FROM linpf.res WHERE prop ='ang' and criteria in ('cor') GROUP BY casename, criteria ) )
SELECT t.* FROM linpf.res as t, tmp 
WHERE t.casename=tmp.casename and t.criteria=tmp.criteria and t.value=tmp.val
ORDER BY t.casename, t.prop, t.criteria;


CREATE TABLE linpf.best_nogw AS
WITH tmp as (
(SELECT casename, 'pf' as prop, criteria, min(value) as val FROM linpf.nogw WHERE prop like 'pf%' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'pf' as prop, criteria, max(value) as val FROM linpf.nogw WHERE prop like 'pf%' and criteria in ('cor') GROUP BY casename, criteria )
UNION
(SELECT casename, 'pt' as prop, criteria, min(value) as val FROM linpf.nogw WHERE prop like 'pt%' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'pt' as prop, criteria, max(value) as val FROM linpf.nogw WHERE prop like 'pt%' and criteria in ('cor') GROUP BY casename, criteria )
UNION
(SELECT casename, 'qf' as prop, criteria, min(value) as val FROM linpf.nogw WHERE prop like 'qf%' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'qf' as prop, criteria, max(value) as val FROM linpf.nogw WHERE prop like 'qf%' and criteria in ('cor') GROUP BY casename, criteria )
UNION
(SELECT casename, 'qt' as prop, criteria, min(value) as val FROM linpf.nogw WHERE prop like 'qt%' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'qt' as prop, criteria, max(value) as val FROM linpf.nogw WHERE prop like 'qt%' and criteria in ('cor') GROUP BY casename, criteria )
UNION
(SELECT casename, 'volt' as prop, criteria, min(value) as val FROM linpf.nogw WHERE prop ='volt' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'volt' as prop, criteria, max(value) as val FROM linpf.nogw WHERE prop ='volt' and criteria in ('cor') GROUP BY casename, criteria )
UNION
(SELECT casename, 'ang' as prop, criteria, min(value) as val FROM linpf.nogw WHERE prop ='ang' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'ang' as prop, criteria, max(value) as val FROM linpf.nogw WHERE prop ='ang' and criteria in ('cor') GROUP BY casename, criteria ) )
SELECT t.* FROM linpf.nogw as t, tmp 
WHERE t.casename=tmp.casename and t.criteria=tmp.criteria and t.value=tmp.val
ORDER BY t.casename, t.prop, t.criteria;


CREATE TABLE linpf.best0 AS
WITH tmp as (
(SELECT casename, 'pf' as prop, criteria, min(value) as val FROM linpf.res WHERE id like '0%' and prop like 'pf%' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'pf' as prop, criteria, max(value) as val FROM linpf.res WHERE id like '0%' and prop like 'pf%' and criteria in ('cor') GROUP BY casename, criteria )
UNION
(SELECT casename, 'pt' as prop, criteria, min(value) as val FROM linpf.res WHERE id like '0%' and prop like 'pt%' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'pt' as prop, criteria, max(value) as val FROM linpf.res WHERE id like '0%' and prop like 'pt%' and criteria in ('cor') GROUP BY casename, criteria )
UNION
(SELECT casename, 'qf' as prop, criteria, min(value) as val FROM linpf.res WHERE id like '0%' and prop like 'qf%' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'qf' as prop, criteria, max(value) as val FROM linpf.res WHERE id like '0%' and prop like 'qf%' and criteria in ('cor') GROUP BY casename, criteria )
UNION
(SELECT casename, 'qt' as prop, criteria, min(value) as val FROM linpf.res WHERE id like '0%' and prop like 'qt%' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'qt' as prop, criteria, max(value) as val FROM linpf.res WHERE id like '0%' and prop like 'qt%' and criteria in ('cor') GROUP BY casename, criteria )
UNION
(SELECT casename, 'volt' as prop, criteria, min(value) as val FROM linpf.res WHERE id like '0%' and prop ='volt' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'volt' as prop, criteria, max(value) as val FROM linpf.res WHERE id like '0%' and prop ='volt' and criteria in ('cor') GROUP BY casename, criteria )
UNION
(SELECT casename, 'ang' as prop, criteria, min(value) as val FROM linpf.res  WHERE id like '0%' and prop ='ang' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'ang' as prop, criteria, max(value) as val FROM linpf.res  WHERE id like '0%' and prop ='ang' and criteria in ('cor') GROUP BY casename, criteria ) )
SELECT t.* FROM linpf.res as t, tmp 
WHERE t.casename=tmp.casename and t.criteria=tmp.criteria and t.value=tmp.val
ORDER BY t.casename, t.prop, t.criteria;

CREATE TABLE linpf.bestnogw0 AS
WITH tmp as (
(SELECT casename, 'pf' as prop, criteria, min(value) as val FROM linpf.nogw WHERE id like '0%' and prop like 'pf%' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'pf' as prop, criteria, max(value) as val FROM linpf.nogw WHERE id like '0%' and prop like 'pf%' and criteria in ('cor') GROUP BY casename, criteria )
UNION
(SELECT casename, 'pt' as prop, criteria, min(value) as val FROM linpf.nogw WHERE id like '0%' and prop like 'pt%' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'pt' as prop, criteria, max(value) as val FROM linpf.nogw WHERE id like '0%' and prop like 'pt%' and criteria in ('cor') GROUP BY casename, criteria )
UNION
(SELECT casename, 'qf' as prop, criteria, min(value) as val FROM linpf.nogw WHERE id like '0%' and prop like 'qf%' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'qf' as prop, criteria, max(value) as val FROM linpf.nogw WHERE id like '0%' and prop like 'qf%' and criteria in ('cor') GROUP BY casename, criteria )
UNION
(SELECT casename, 'qt' as prop, criteria, min(value) as val FROM linpf.nogw WHERE id like '0%' and prop like 'qt%' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'qt' as prop, criteria, max(value) as val FROM linpf.nogw WHERE id like '0%' and prop like 'qt%' and criteria in ('cor') GROUP BY casename, criteria )
UNION
(SELECT casename, 'volt' as prop, criteria, min(value) as val FROM linpf.nogw WHERE id like '0%' and prop ='volt' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'volt' as prop, criteria, max(value) as val FROM linpf.nogw WHERE id like '0%' and prop ='volt' and criteria in ('cor') GROUP BY casename, criteria )
UNION
(SELECT casename, 'ang' as prop, criteria, min(value) as val FROM linpf.nogw  WHERE id like '0%' and prop ='ang' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'ang' as prop, criteria, max(value) as val FROM linpf.nogw  WHERE id like '0%' and prop ='ang' and criteria in ('cor') GROUP BY casename, criteria ) )
SELECT t.* FROM linpf.nogw as t, tmp 
WHERE t.casename=tmp.casename and t.criteria=tmp.criteria and t.value=tmp.val
ORDER BY t.casename, t.prop, t.criteria;

CREATE TABLE linpf.best1 AS
WITH tmp as (
(SELECT casename, 'pf' as prop, criteria, min(value) as val FROM linpf.res WHERE id like '1%' and prop like 'pf%' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'pf' as prop, criteria, max(value) as val FROM linpf.res WHERE id like '1%' and prop like 'pf%' and criteria in ('cor') GROUP BY casename, criteria )
UNION
(SELECT casename, 'pt' as prop, criteria, min(value) as val FROM linpf.res WHERE id like '1%' and prop like 'pt%' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'pt' as prop, criteria, max(value) as val FROM linpf.res WHERE id like '1%' and prop like 'pt%' and criteria in ('cor') GROUP BY casename, criteria )
UNION
(SELECT casename, 'qf' as prop, criteria, min(value) as val FROM linpf.res WHERE id like '1%' and prop like 'qf%' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'qf' as prop, criteria, max(value) as val FROM linpf.res WHERE id like '1%' and prop like 'qf%' and criteria in ('cor') GROUP BY casename, criteria )
UNION
(SELECT casename, 'qt' as prop, criteria, min(value) as val FROM linpf.res WHERE id like '1%' and prop like 'qt%' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'qt' as prop, criteria, max(value) as val FROM linpf.res WHERE id like '1%' and prop like 'qt%' and criteria in ('cor') GROUP BY casename, criteria )
UNION
(SELECT casename, 'volt' as prop, criteria, min(value) as val FROM linpf.res WHERE id like '1%' and prop ='volt' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'volt' as prop, criteria, max(value) as val FROM linpf.res WHERE id like '1%' and prop ='volt' and criteria in ('cor') GROUP BY casename, criteria )
UNION
(SELECT casename, 'ang' as prop, criteria, min(value) as val FROM linpf.res  WHERE id like '1%' and prop ='ang' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'ang' as prop, criteria, max(value) as val FROM linpf.res  WHERE id like '1%' and prop ='ang' and criteria in ('cor') GROUP BY casename, criteria ) )
SELECT t.* FROM linpf.res as t, tmp 
WHERE t.casename=tmp.casename and t.criteria=tmp.criteria and t.value=tmp.val
ORDER BY t.casename, t.prop, t.criteria;

CREATE TABLE linpf.bestnogw1 AS
WITH tmp as (
(SELECT casename, 'pf' as prop, criteria, min(value) as val FROM linpf.nogw WHERE id like '1%' and prop like 'pf%' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'pf' as prop, criteria, max(value) as val FROM linpf.nogw WHERE id like '1%' and prop like 'pf%' and criteria in ('cor') GROUP BY casename, criteria )
UNION
(SELECT casename, 'pt' as prop, criteria, min(value) as val FROM linpf.nogw WHERE id like '1%' and prop like 'pt%' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'pt' as prop, criteria, max(value) as val FROM linpf.nogw WHERE id like '1%' and prop like 'pt%' and criteria in ('cor') GROUP BY casename, criteria )
UNION
(SELECT casename, 'qf' as prop, criteria, min(value) as val FROM linpf.nogw WHERE id like '1%' and prop like 'qf%' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'qf' as prop, criteria, max(value) as val FROM linpf.nogw WHERE id like '1%' and prop like 'qf%' and criteria in ('cor') GROUP BY casename, criteria )
UNION
(SELECT casename, 'qt' as prop, criteria, min(value) as val FROM linpf.nogw WHERE id like '1%' and prop like 'qt%' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'qt' as prop, criteria, max(value) as val FROM linpf.nogw WHERE id like '1%' and prop like 'qt%' and criteria in ('cor') GROUP BY casename, criteria )
UNION
(SELECT casename, 'volt' as prop, criteria, min(value) as val FROM linpf.nogw WHERE id like '1%' and prop ='volt' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'volt' as prop, criteria, max(value) as val FROM linpf.nogw WHERE id like '1%' and prop ='volt' and criteria in ('cor') GROUP BY casename, criteria )
UNION
(SELECT casename, 'ang' as prop, criteria, min(value) as val FROM linpf.nogw  WHERE id like '1%' and prop ='ang' and criteria in ('del', 'max', 'avg') GROUP BY casename, criteria )
UNION 
(SELECT casename, 'ang' as prop, criteria, max(value) as val FROM linpf.nogw  WHERE id like '1%' and prop ='ang' and criteria in ('cor') GROUP BY casename, criteria ) )
SELECT t.* FROM linpf.nogw as t, tmp 
WHERE t.casename=tmp.casename and t.criteria=tmp.criteria and t.value=tmp.val
ORDER BY t.casename, t.prop, t.criteria;


-- ====================================================================
-- From new formulation (tables are linpf.vt, linpf.flows, linpf.cases)
-- ====================================================================

CREATE VIEW linpf.vtmin AS
WITH tmp as (
	(SELECT caseid, prop, criteria, min(round(value::NUMERIC, 6)) as minv FROM linpf.vt WHERE criteria in ('del', 'max', 'avg', 'rms') GROUP BY caseid, prop, criteria)
	UNION
	(SELECT caseid, prop, criteria, max(round(value::NUMERIC,6)) as minv FROM linpf.vt WHERE criteria in ('cor') GROUP BY caseid, prop, criteria) )
SELECT * from tmp ORDER BY caseid, prop, criteria;

CREATE MATERIALIZED VIEW linpf.vtdcount AS
WITH tmp AS (
	SELECT a.*, b.minv FROM linpf.vt a, linpf.vtmin b WHERE a.caseid=b.caseid and a.prop=b.prop and a.criteria=b.criteria),
tmp2 AS (SELECT distinct id FROM linpf.vt)
SELECT tmp2.id , (SELECT count(*) from tmp WHERE round(tmp.value::NUMERIC,6)=tmp.minv and tmp.id=tmp2.id) 
FROM tmp2 ORDER BY id;

--CREATE VIEW linpf.vtidmin AS
--SELECT a.* FROM linpf.vt a, linpf.vtmin b WHERE a.caseid=b.caseid and a.prop=b.prop and a.criteria=b.criteria and round(a.value::NUMERIC, 6)=b.minv
--ORDER BY a.caseid, a.prop, a.criteria, a.id;

CREATE VIEW linpf.vt_casesum AS
SELECT id, prop, criteria, round(sum(value)::NUMERIC, 6) vsum FROM linpf.vt GROUP BY id, prop, criteria ORDER BY prop, criteria, vsum;

CREATE VIEW linpf.vt_minsum AS
(SELECT  a.prop, a.criteria, min(vsum) vmin FROM linpf.vt_casesum a WHERE a.criteria in ('del', 'max', 'avg', 'rms') GROUP BY a.prop, a.criteria ORDER BY a.prop, a.criteria)
UNION
(SELECT  b.prop, b.criteria, max(vsum) vmin FROM linpf.vt_casesum b WHERE b.criteria='cor' GROUP BY b.prop, b.criteria ORDER BY b.prop, b.criteria);

CREATE VIEW linpf.vtbestid AS
SELECT t.* FROM  linpf.vt_casesum t, linpf.vt_minsum t2 WHERE t.prop=t2.prop and t.criteria=t2.criteria and t.vsum=t2.vmin ORDER BY t.prop, t.criteria;

CREATE VIEW linpf.vt_corsum AS
SELECT id, sum(value) cv FROM linpf.vt WHERE criteria='cor' GROUP BY id;

CREATE VIEW linpf.vt_noncorsum AS
SELECT id, sum(value) FROM linpf.vt WHERE criteria in ('del', 'max', 'avg', 'rms') GROUP BY id;

CREATE VIEW linpf.vt_idsum AS
SELECT a.id, a.sum - b.cv as v FROM linpf.vt_noncorsum a, linpf.vt_corsum b WHERE a.id=b.id;

CREATE VIEW linpf.vt_noncorsum2 AS
SELECT id, sum(value) FROM linpf.vt WHERE criteria in ('max', 'avg', 'rms') GROUP BY id;

CREATE VIEW linpf.vt_idsum2 AS
SELECT a.id, a.sum - b.cv as v FROM linpf.vt_noncorsum2 a, linpf.vt_corsum b WHERE a.id=b.id;

CREATE MATERIALIZED VIEW linpf.idmap AS
	
