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
--Snapshot to be analized
DECLARE @snapnum INT
--Hubble constant
DECLARE @H FLOAT

--Asignment of variables
SET @x = 0
SET @y = 0
SET @z = 0
SET @id = 1234
SELECT @hx = M.x,
       @hy = M.y,
       @hz = M.z,
       @vx = M.velX,
       @vy = M.velY,
       @vz = M.velZ
FROM MPAHalo M
WHERE M.haloID = @id

--Selecciona el primer neighbour del halo

WHERE D.
