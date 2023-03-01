from astropy.table import QTable


class Quality_Control(object):

	def __init__(self, file_name):

		self._harvest_tables(file_name)

	def _harvest_tables(self, file_name):

		with open(file_name, 'r') as f:
			d, header = self._harvest_data(f)
			f.close()

		table_dict = {key: {} for key in d}
		for key in d:
			for i, colname in enumerate(header):
				table_dict[key][colname] = [row[i] for row in d[key]]

		for key in table_dict:
			t = QTable([table_dict[key][k] for k in table_dict[key]], names=[k for k in table_dict[key]])
			setattr(self, "%s_table"%key, t)
	
	def _harvest_data(self, f):

		d = {"object": [],"flat": [],"arc": [],}

		i = 0
		in_table = False ; current_table = None ; header = None
		for line in f:
			if "IMAGE" in line.split():
				header = line.split()
				in_table = True
				current_table = list(d.keys())[0]
			elif in_table and line[0].isalpha() and header is not None:
				in_table = False
			elif not in_table and line[0]=='-' and header is not None:
				i += 1
				if i > 2: break
				current_table = list(d.keys())[i]
				in_table=True
			elif in_table and len(line.split())>4:
				d[current_table] += [line.split()]

		return d, header


