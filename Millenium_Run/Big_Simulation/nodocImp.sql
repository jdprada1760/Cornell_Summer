--Position of the point of view
DECLARE @x FLOAT
DECLARE @y FLOAT
DECLARE @z FLOAT
DECLARE @H FLOAT
DECLARE @lh float
DECLARE @neigh INT
DECLARE @limxh INT
DECLARE @limyh INT
DECLARE @limzh INT
DECLARE @limMxh INT
DECLARE @limMyh INT
DECLARE @limMzh INT
DECLARE @snapnum INT
DECLARE @nb INT
SET @nb = 3
SET @snapnum = 63
SET @limxh = %x
SET @limyh = %y
SET @limzh = %z
SET @limMxh = $x
SET @limMyh = $y
SET @limMzh = $z
SET @lh = 0.7
SET @snapnum = 63
SET @neigh = 3
SET @x = 0
SET @y = 0
SET @z = 0
SET @H = 100
SELECT *
FROM    (
        SELECT *, ROW_NUMBER() OVER ( PARTITION BY f.haloID ORDER BY f.SIGMA ASC ) AS RN
        FROM  (SELECT   ((( POWER(M.x - @x,2) + POWER(M.y - @y,2) + POWER(M.z - @z,2) )
                        *( POWER(D.x - @x,2) + POWER(D.y - @y,2) + POWER(D.z - @z,2) )
                        /POWER((M.x - @x)*(D.x - @x) + (M.y - @y)*(D.y - @y) + (M.z - @z)*(D.z - @z), 2))
                        -1)*( POWER(M.x - @x,2) + POWER(M.y - @y,2) + POWER(M.z - @z,2)) AS SIGMA,
                         M.np AS NP,
                         M.haloID AS haloID
              FROM (SELECT galax.*, MaH.x, MaH.y, MaH.z, MaH.velX, MaH.velY, MaH.velZ, MaH.np
                    FROM  (SELECT Ma.haloID, Ma.galaxyID
                          FROM MPAGalaxies..DeLucia2006a Ma
                          WHERE Ma.snapNum = @snapnum
                                AND Ma.ix >= FLOOR(@limxh/10) AND Ma.iy >= FLOOR(@limyh/10) AND Ma.iz >= FLOOR(@limzh/10)
                                AND Ma.ix < FLOOR(@limMxh/10) AND Ma.iy < FLOOR(@limMyh/10) AND Ma.iz < FLOOR(@limMzh/10)
                                AND Ma.type = 0) galax
                    INNER JOIN (SELECT halo.haloID, halo.x, halo.y, halo.z, halo.velX, halo.velY, halo.velZ, halo.np
                               FROM MPAHaloTrees..MHalo halo
                               WHERE halo.snapNum = @snapnum
                                      AND halo.x >= @limxh - 5 AND halo.y >= @limyh - 5 AND halo.z >= @limzh - 5
                                      AND halo.x < @limMxh + 5 AND halo.y < @limMyh + 5 AND halo.z < @limMzh + 5 ) MaH
                    ON MaH.haloID = galax.haloID ) M,
                   (SELECT a.*
                     FROM (SELECT D.galaxyID, D.ix, D.iy, D.iz, D.x, D.y, D.z, D.velX, D.velY, D.velZ
                           FROM MPAGalaxies..DeLucia2006a D
                           WHERE  D.snapNum = @snapnum
                                 AND D.ix >= FLOOR(@limxh/10) - @nb AND D.iy >= FLOOR(@limyh/10) - @nb AND D.iz >= FLOOR(@limzh/10) - @nb
                                 AND D.ix <= FLOOR(@limMxh/10) + @nb AND D.iy <= FLOOR(@limMyh/10) + @nb AND D.iz <= FLOOR(@limMzh/10) + @nb
                           ) a
                     INNER JOIN (SELECT L.galaxyID
                                FROM MPAGalaxies..Delucia2006a_SDSS2MASS L
                                WHERE L.snapNum = @snapnum AND L.r_sdss < -19) S
                     ON S.galaxyID = a.galaxyID ) D
              WHERE   D.galaxyID != M.galaxyID AND
                      ABS(((@H*(D.x - M.x)/@lh + D.velX - M.velX)*(M.x - @x)+
                         (@H*(D.y - M.y)/@lh + D.velY - M.velY)*(M.y - @y) +
                         (@H*(D.z - M.z)/@lh + D.velZ - M.velZ)*(M.z - @z))/SQRT( POWER(M.x - @x,2) + POWER(M.y - @y,2) + POWER(M.z - @z,2) )
                        )  < 500) f
      ) g
WHERE RN = @neigh
