--This query has into account haloes with not enough neighbouring galaxies to have nth closest neighbor
--Position of the point of view
DECLARE @x FLOAT
DECLARE @y FLOAT
DECLARE @z FLOAT
--Snapshot to be analyzed
DECLARE @snapnum INT
--Hubble constant
DECLARE @H FLOAT
-- Little h
DECLARE @lh float
--Neighbour
DECLARE @neigh INT
--Limits the box of halos (this time uses the number associated to each position pointer)
DECLARE @limx INT
DECLARE @limy INT
DECLARE @limz INT
DECLARE @limMx INT
DECLARE @limMy INT
DECLARE @limMz INT
--Limits the size of the cube of galaxies for each halo
DECLARE @nb INT


--Asignment of variables
SET @nb = $n
SET @lh = 0.7
SET @snapnum = 63
SET @neigh = 4
SET @x = 0
SET @y = 0
SET @z = 0
SET @limx = %x
SET @limy = %y
SET @limz = %z
SET @limMx = $x
SET @limMy = $y
SET @limMz = $z
SET @H = 100


SELECT  POWER(@nb*10,2) AS SIGMA, MAX(f.NP) AS NP, g.haloID AS haloID, COUNT(f.haloID) AS RN
FROM
    --Gets the haloes and galaxies IDs that are near (nb boxes close)
    (SELECT
               M.np AS NP,
               M.haloID AS haloID
    FROM MPAHalo M,
         DeLucia2006a D
    WHERE --Selects snapshot
          D.snapnum = @snapnum
          AND M.snapNum = @snapnum
          --Omits galaxies with r abs magnitude < -19
          AND D.mag_r < -19
          --Restricts the box of haloes
          AND M.ix >= @limx AND M.iy >= @limy AND M.iz >= @limz
          AND M.ix <= @limMx AND M.iy <= @limMy AND M.iz <= @limMz
          --Restricts the box of galaxies for each halo
          AND D.ix >= M.ix - @nb AND D.iy >= M.iy - @nb AND D.iz >= M.iz - @nb
          AND D.ix <= M.ix + @nb AND D.iy <= M.iy + @nb AND D.iz <= M.iz + @nb
          --Selects only galaxies that fulfill the velocity cut
          AND ABS(
               ((@H*(D.x - M.x)/@lh + D.velX - M.velX)*(M.x - @x)+
               (@H*(D.y - M.y)/@lh + D.velY - M.velY)*(M.y - @y) +
               (@H*(D.z - M.z)/@lh + D.velZ - M.velZ)*(M.z - @z))/SQRT( POWER(M.x - @x,2) + POWER(M.y - @y,2) + POWER(M.z - @z,2) )

                 )  < 500
) f
RIGHT JOIN (SELECT M.haloID AS haloID
            FROM MPAHalo M
            WHERE M.snapNum = @snapnum
            AND M.ix >= @limx AND M.iy >= @limy AND M.iz >= @limz
            AND M.ix <= @limMx AND M.iy <= @limMy AND M.iz <= @limMz) g
ON f.haloID = g.haloID
GROUP BY g.haloID
ORDER BY COUNT(f.haloID)
HAVING COUNT(f.haloID) < @neigh
