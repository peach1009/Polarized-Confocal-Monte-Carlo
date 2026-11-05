function parsaveLastPhotonData(filename, lastPos, lastScat, lastStokes)
    save(filename, 'lastPos', 'lastScat', 'lastStokes', '-v7.3');
end