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
--Minimum and Maximum masses
DECLARE @minM FLOAT
DECLARE @maxM FLOAT
--Number of bins in the histogram
DECLARE @nbins INT

--Assignment of values
SET @snapnum = 63
SET @nbins = 10
SET @l = 0  
SET @bsize = 10  
SELECT @posx = @l*RAND()  
SELECT @posy = @l*RAND()  
SELECT @posz = @l*RAND()  

--Selects the minimum and maximum masses
SELECT @minM = MIN(D.np),
       @maxM = MAX(D.np)
FROM MPAHalo D
WHERE D.x > @posx AND D.x < @posx + @bsize  
  AND D.y > @posx AND D.y < @posy + @bsize  
  AND D.z > @posx AND D.z < @posz + @bsize
  AND D.snapnum = @snapnum

--Defines the interval of the histogram
DECLARE @interv FLOAT
SET @interv = ( @maxM - @minM ) / @nbins

--Selects the histogram
SELECT @interv*(FLOOR(D.np/@interv)) AS mag,
       COUNT(*) AS NUM
FROM MPAHalo D
WHERE D.x > @posx AND D.x < @posx + @bsize  
  AND D.y > @posx AND D.y < @posy + @bsize  
  AND D.z > @posx AND D.z < @posz + @bsize
  AND D.snapnum = @snapnum
GROUP BY @interv*(FLOOR(D.np/@interv))
ORDER BY mag