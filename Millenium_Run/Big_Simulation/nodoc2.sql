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
SET @snapnum = 63
SET @limxh = %x
SET @limyh = %y
SET @limzh = %z
SET @limMxh = $x
SET @limMyh = $y
SET @limMzh = $z
SET @lh = 0.7
SET @snapnum = 63
SET @neigh = 1
SET @x = 0
SET @y = 0
SET @z = 0
SET @H = 100
SELECT *
FROM    (
        SELECT *, ROW_NUMBER() OVER ( PARTITION BY f.haloID ORDER BY f.SIGMA ASC ) AS RN
        FROM  (SELECT (
                         POWER(((D.x - @x)*(M.x - @x)
                                + (D.y - @y)*(M.y - @y)
                                + (D.z - @z)*(M.z - @z)
                                )*(M.x - @x)/( POWER(M.x - @x,2) + POWER(M.y - @y,2) + POWER(M.z - @z,2) )
                                - (D.x - @x), 2 )+
                         POWER(((D.x - @x)*(M.x - @x)
                                + (D.y - @y)*(M.y - @y)
                                + (D.z - @z)*(M.z - @z))*(M.y - @y)/( POWER(M.x - @x,2) + POWER(M.y - @y,2) + POWER(M.z - @z,2) )
                                - (D.y - @y)  , 2 )+
                         POWER(((D.x - @x)*(M.x - @x)
                                + (D.y - @y)*(M.y - @y)
                                + (D.z - @z)*(M.z - @z))*(M.z - @z)/( POWER(M.x - @x,2) + POWER(M.y - @y,2) + POWER(M.z - @z,2) )
                                - (D.z - @z), 2 )) AS SIGMA,
                         M.np AS NP,
                         M.haloID AS haloID
              FROM (SELECT Ma.haloID, Ma.x, Ma.y, Ma.z, Ma.velX, Ma.velY, Ma.velZ, Ma.np
                    FROM MPAHaloTrees..MHalo Ma
                    WHERE Ma.snapNum = @snapnum
                          AND Ma.x >= @limxh AND Ma.y >= @limyh AND Ma.z >= @limzh
                          AND Ma.x < @limMxh AND Ma.y < @limMyh AND Ma.z < @limMzh) M,
                   (SELECT a.*
                     FROM (SELECT D.galaxyID, D.ix, D.iy, D.iz, D.x, D.y, D.z
                           FROM MPAGalaxies..DeLucia2006a D
                           WHERE  D.snapNum = @snapnum) a
                     RIGHT JOIN (SELECT L.galaxyID
                                FROM MPAGalaxies..Delucia2006a_SDSS2MASS L
                                WHERE L.snapNum = @snapnum AND L.K_2mass < -25) S
                     ON S.galaxyID = a.galaxyID ) D
      )f
      ) g

WHERE RN = @neigh
