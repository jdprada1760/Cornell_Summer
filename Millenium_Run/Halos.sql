--Position of a corner of the box
DECLARE @posx FLOAT
DECLARE @posy FLOAT
DECLARE @posz FLOAT
--Scale to generate random number
DECLARE @l FLOAT
--Size of the box
DECLARE @bsize FLOAT

SET @l = 100
SET @bsize = 10
SELECT @posx = @l*RAND()
SELECT @posy = @l*RAND()
SELECT @posz = @l*RAND()

SELECT D.haloId
FROM MPAHalo D
WHERE D.x > @posx AND D.x < @posx + @bsize  
  AND D.y > @posx AND D.y < @posy + @bsize  
  AND D.z > @posx AND D.z < @posz + @bsize 
;
