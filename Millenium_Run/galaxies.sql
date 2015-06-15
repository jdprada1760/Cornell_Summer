--Position of a corner of the box
DECLARE @posx FLOAT
DECLARE @posy FLOAT
DECLARE @posz FLOAT
--Snapshot to be analized
DECLARE @snapnum INT
--Scale to generate random number
DECLARE @l FLOAT
--Size of the box
DECLARE @bsize FLOAT
--Minimum and Maximum logarithms of Masses and Luminosities
DECLARE @minM FLOAT
DECLARE @maxM FLOAT
DECLARE @minL FLOAT
DECLARE @maxL FLOAT
--Number of bins in the histogram
DECLARE @nbins INT
--Mass (in solar masses) unit of virial mass
DECLARE @mp FLOAT
--Solar luminosity (erg/s)
DECLARE @lsun FLOAT
--Solar visual absolute magnitude
DECLARE @sunmv FLOAT

--Assignment of values
SET @sunmv = 4.78
SET @lsun = 3.846E+33
SET @mp = 1E+10
SET @snapnum = 63
SET @nbins = 20
SET @l = 0
SET @bsize = 10000000
SELECT @posx = @l*RAND()
SELECT @posy = @l*RAND()
SELECT @posz = @l*RAND()

--Selects the minimum and maximum masses
SELECT @minM = MIN(LOG10(@mp*D.mvir)),
       @maxM = MAX(LOG10(@mp*D.mvir))
FROM DeLucia2006a D
WHERE D.x > @posx AND D.x < @posx + @bsize
  AND D.y > @posx AND D.y < @posy + @bsize
  AND D.z > @posx AND D.z < @posz + @bsize
  AND D.snapnum = @snapnum

--Selects the minimum and maximum magnitudes
SELECT @minL = MIN((@sunmv - D.mag_v)/2.5),
       @maxL = MAX((@sunmv - D.mag_v)/2.5)
FROM DeLucia2006a D
WHERE D.x > @posx AND D.x < @posx + @bsize
  AND D.y > @posx AND D.y < @posy + @bsize
  AND D.z > @posx AND D.z < @posz + @bsize
  AND D.snapnum = @snapnum
  AND D.mag_v != 99.0

--Defines the interval of the histogram for masses and luminosities
DECLARE @intervm FLOAT
DECLARE @intervL FLOAT
SET @intervm = ( @maxM - @minM ) / @nbins
SET @intervL = ( @maxL - @minL ) / @nbins

--Selects the histogram logarithmic scale for masses
SELECT @intervm*(FLOOR((LOG10(D.mvir*@mp)-@minM)/@intervm))+@minM AS logM,
       COUNT(*) AS NUM_M,
       @minM AS MIN_M,
       @maxM AS MAX_M,
       @nbins AS NBINS
FROM DeLucia2006a D
WHERE D.x > @posx AND D.x < @posx + @bsize
      AND D.y > @posx AND D.y < @posy + @bsize
      AND D.z > @posx AND D.z < @posz + @bsize
      AND D.snapnum = @snapnum
GROUP BY @intervm*(FLOOR((LOG10(D.mvir*@mp)-@minM)/@intervm))+@minM
ORDER BY logM



--Selects the histogram logarithmic scale for luminosities
SELECT @intervL*(FLOOR((((@sunmv - D.mag_v)/2.5)-@minL)/@intervL))+@minL AS logL,
       COUNT(*) AS NUM_L,
       @minL AS MIN_L,
       @maxL AS MAX_L,
       @nbins AS NBINS
FROM DeLucia2006a D
WHERE D.x > @posx AND D.x < @posx + @bsize
      AND D.y > @posx AND D.y < @posy + @bsize
      AND D.z > @posx AND D.z < @posz + @bsize
      AND D.snapnum = @snapnum
      AND D.mag_v != 99.0
GROUP BY @intervL*(FLOOR((((@sunmv - D.mag_v)/2.5)-@minL)/@intervL))+@minL
ORDER BY logL
