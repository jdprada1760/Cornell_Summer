--Position of a corner of the box  
DECLARE @posx FLOAT  
DECLARE @posy FLOAT  
DECLARE @posz FLOAT  
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

--Prints the minimum and maximum masses
/*SELECT 'Maximum Mass: ' + CAST( @maxM AS VARCHAR(15) ) + CHAR(13) + ' Minimum Mass: ' + CAST( @minM AS VARCHAR(15) ) + CHAR(13) + CHAR(13) + CHAR(13)'*/

--Creates the histogram table
DECLARE @histo TABLE(
  Lmass FLOAT,
  Rmass FLOAT,
  nHaloes INT
);

--Fills the table with the values
DECLARE @i INT
DECLARE @cnt INT
DECLARE @lmass FLOAT
DECLARE @rmass FLOAT
DECLARE @interv FLOAT
SET @interv = ( @maxM - @minM ) / @nbins
SET @lmass = @minM
SET @rmass = @lmass + @interv
SET @i = @nbins
WHILE ( @i > 0 )
BEGIN
--Counts the haloes within the mass interval given by rmass and lmass
SELECT @cnt = COUNT(D.haloId)
       FROM MPAHalo D
       WHERE D.x > @posx AND D.x < @posx + @bsize  
         AND D.y > @posx AND D.y < @posy + @bsize  
	 AND D.z > @posx AND D.z < @posz + @bsize
	 AND D.np >= @lmass
	 AND D.np < @rmass

--Insert the values into the histo table
INSERT INTO @histo( Lmass, Rmass, nHaloes )
      VALUES( @lmass, @rmass, @cnt )

--Update the index and the interval limits
SET @lmass = @rmass
SET @rmass = @rmass + @interv
SET @i = @i - 1

END

SELECT H.nHaloes, H.Lmass, H.Rmass
FROM @histo H
;