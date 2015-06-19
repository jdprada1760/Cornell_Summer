--Position of the point of view
DECLARE @x FLOAT
DECLARE @y FLOAT
DECLARE @z FLOAT
--Halo ID, position and peculiar velocity
DECLARE @id INT
DECLARE @hx FLOAT
DECLARE @hy FLOAT
DECLARE @hz FLOAT
DECLARE @vx FLOAT
DECLARE @vy FLOAT
DECLARE @vz FLOAT
--Normalized vector from pov to halo
DECLARE @dx FLOAT
DECLARE @dy FLOAT
DECLARE @dz FLOAT
--Snapshot to be analized
DECLARE @snapnum INT
--Hubble constant
DECLARE @H FLOAT
--Neighbour
DECLARE @neigh INT

--Asignment of variables
SET @neigh = 4
SET @x = 250
SET @y = 250
SET @z = 250
SET @H = 100
SET @id = 1234
SELECT @hx = M.x,
       @hy = M.y,
       @hz = M.z,
       @vx = M.velX,
       @vy = M.velY,
       @vz = M.velZ
FROM MPAHalo M
WHERE M.haloID = @id
SET @dx = (@hx - @x)/SQRT( POWER(@hx - @x,2) + POWER(@hy - @y,2) + POWER(@hz - @z,2) )
SET @dy = (@hy - @y)/SQRT( POWER(@hx - @x,2) + POWER(@hy - @y,2) + POWER(@hz - @z,2) )
SET @dz = (@hz - @z)/SQRT( POWER(@hx - @x,2) + POWER(@hy - @y,2) + POWER(@hz - @z,2) )




--Makes an ordered list distances to galaxies
SELECT *, SQRT(POWER(b.vx,2) + POWER(b.vy,2)+ POWER(b.vz,2))
FROM
	(SELECT TOP 10 ROW_NUMBER() OVER (ORDER BY a.DIST) AS NUM,
       	a.DIST,
	a.vx,
	a.vy,
	a.vz
	FROM
		(SELECT SQRT(
		       POWER(  ((D.x - @x)*@dx + (D.y - @y)*@dy + (D.z - @z)*@dz)*@dx - (D.x - @x)  , 2 )  +
		       POWER(  ((D.x - @x)*@dx + (D.y - @y)*@dy + (D.z - @z)*@dz)*@dy - (D.y - @y)  , 2 )  +
	      	       POWER(  ((D.x - @x)*@dx + (D.y - @y)*@dy + (D.z - @z)*@dz)*@dz - (D.z - @z)  , 2 )  
	                    ) AS DIST,
		       D.velX AS vx,
		       D.velY AS vy,
		       D.velZ AS vz
	 	FROM DeLucia2006a D
	 	WHERE --Omits galaxies with r abs magnitude < -19
	       	D.mag_r < -19
	      	 --Only counts galxies with relative velocities to the halo < 500 Km/s
	       	 AND ABS(  
	  	       --Radial component of the total velocity of the galaxy
	       	       @dx*(@H*(D.x-@x) + D.velX) + @dy*(@H*(D.y-@y) + D.velY) + @dz*(@H*(D.z-@z) + D.velZ)
		       --Radial component of the total velocity of the halo
		       -(@dx*(@H*(@hx-@x) + @vx) + @dy*(@H*(@hy-@y) + @vy) + @dz*(@H*(@hz-@z) + @vz))

		      ) < 500
		) a
	ORDER BY DIST
	)b
WHERE b.NUM = @neigh