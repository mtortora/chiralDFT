from VolumeViewer import volume


drlist = volume.volume_list()

for dr in drlist: dr.remove_surfaces()