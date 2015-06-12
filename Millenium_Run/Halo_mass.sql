--Position of a corner of the box
DECLARE @posx FLOAT
DECLARE @posy FLOAT
DECLARE @posz FLOAT
--Scale to generate random number
DECLARE @l FLOAT
--Size of the box
DECLARE @bsize FLOAT
--Snapshot to be analized
DECLARE @snapnum INT

--Assignment of values
SET @snapnum = 63
SET @l = 0
SET @bsize = 10
SELECT @posx = @l*RAND()
SELECT @posy = @l*RAND()
SELECT @posz = @l*RAND()

--Prints the minimum and maximum masses
/*SELECT 'Maximum Mass: ' + CAST( @maxM AS VARCHAR(15) ) + CHAR(13) + ' Minimum Mass: ' + CAST( @minM AS VARCHAR(15) ) + CHAR(13) + CHAR(13) + CHAR(13)'*/

--Gets the table of the halo masses
SELECT D.haloId, D.np
FROM MPAHalo D
WHERE D.x > @posx AND D.x < @posx + @bsize
      AND D.y > @posx AND D.y < @posy + @bsize
      AND D.z > @posx AND D.z < @posz + @bsize
      AND D.snapnum = @snapnum
