from ..athena.athena_read import athinput
from numpy import logspace, log10, zeros
from netCDF4 import Dataset
from .utils import get_rt_bands, get_ray_out
import re, subprocess

def create_inputs(tmpfile, args):
    with open(tmpfile, 'r') as file:
      tmpinp = file.read()

    pmax, pmin, np = tuple(map(float, args['plevel'].split(':')))
    np = int(np)
    plevel = logspace(log10(pmax), log10(pmin), np)
    plevel = ['%10.2f'%x for x in plevel]

    if args['tem'] == '0':
        Tp = ['0.']*np
    else:
        Tp = args['tem'].split(' ')
    assert(len(Tp) == np)

    if args['nh3'] == '0':
        NH3p  = ['0.']*np
    else:
        NH3p = args['nh3'].split(' ')
    assert(len(NH3p) == np)

  # adjust minimum number of walkers
    #args['nwalker'] = str(max(int(args['nwalker']), 2*np))

    var = [x for x in args['var']]

    name = tmpfile.split('.')[0]
    if args['output'] != '':
        name += '-' + args['output']

    inpfile = re.sub('\[problem_id\]', name, tmpinp)
    inpfile = re.sub('\[logname\]', name, inpfile)
    inpfile = re.sub('\[obsname\]', args['obs'], inpfile)
    inpfile = re.sub('\[T0\]', args['T0'], inpfile)
    inpfile = re.sub('\[Tmin\]', args['Tmin'], inpfile)
    inpfile = re.sub('\[qH2O\]', args['qH2O'], inpfile)
    inpfile = re.sub('\[qNH3\]', args['qNH3'], inpfile)
    inpfile = re.sub('\[Tstd\]', args['sT'], inpfile)
    inpfile = re.sub('\[Xstd\]', args['sNH3'], inpfile)
    inpfile = re.sub('\[Tlen\]', args['zT'], inpfile)
    inpfile = re.sub('\[Xlen\]', args['zNH3'], inpfile)
    inpfile = re.sub('\[plevel\]', ' '.join(plevel), inpfile)
    inpfile = re.sub('\[pmin\]', args['pmin'], inpfile)
    inpfile = re.sub('\[pmax\]', args['pmax'], inpfile)
    inpfile = re.sub('\[nwalker\]', args['nwalker'], inpfile)
    inpfile = re.sub('\[nlim\]', args['nlim'], inpfile)
    inpfile = re.sub('\[variables\]', ' '.join(var), inpfile)
    inpfile = re.sub('\[nodes\]', str(4*int(args['nodes'])), inpfile)
    inpfile = re.sub('\[Tp\]', ' '.join(Tp), inpfile)
    inpfile = re.sub('\[NH3p\]', ' '.join(NH3p), inpfile)
    inpfile = re.sub('\[grav\]', '-' + args['grav'], inpfile)
    inpfile = re.sub('\[lat\]', args['lat'], inpfile)
    inpfile = re.sub('\[rgradt\]', args['rgradt'], inpfile)
    inpfile = re.sub('\[metallicity\]', args['metallicity'], inpfile)
    inpfile = re.sub('\[karpowicz_scale\]', args['karpowicz_scale'], inpfile)
    inpfile = re.sub('\[hanley_power\]', args['hanley_power'], inpfile)
    if args['d']:
        inpfile = re.sub('\[diff\]', 'true', inpfile)
    else:
        inpfile = re.sub('\[diff\]', 'false', inpfile)

    with open(name + '.inp', 'w') as file:
        file.write(inpfile)
    print('Input file written to %s.inp\n' % name)
    return name + '.inp'

def run_forward(exefile, inpfile):
    script = ['./' + exefile, '-i', inpfile]
    process = subprocess.Popen(script,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE)
    while True:
        output = process.stdout.readline()
        #err = process.stderr.readline()
        #if (err != ''):
        #  raise Exception(err.decode('UTF-8'))
        if output == b'' and process.poll() is not None:
            break
        if output:
            print(output.decode('UTF-8'), end = '')
    process.poll()
    #print(out.decode('UTF-8'), end = '\r')
    #print(err.decode('UTF-8'))

    out, err = subprocess.Popen('./combine.py',
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE).communicate()
    print(out.decode('UTF-8'), end = '')
    print(err.decode('UTF-8'), end = '')

    inp = athinput(inpfile)
    return inp['job']['problem_id']

def write_observation(inpfile, datafile):
    freq = get_rt_bands(inpfile)[:,0]
    num_bands = len(freq)
    amu, aphi = get_ray_out(inpfile)
    num_dirs = len(amu)

# read radiation toa
#datafile = 'vla_ideal_js-main.nc'
    data = Dataset(datafile, 'r')
    tb = zeros((4, num_bands, num_dirs))
    for i in range(num_bands):
        for j in range(num_dirs):
            tb[0,i,j] = data['b%dtoa%d' % (i+1,j+1)][0,0,0]
            tb[1,i,j] = data['b%dtoa%d' % (i+1,j+1)][0,1,0]
            tb[2,i,j] = data['b%dtoa%d' % (i+1,j+1)][0,2,0]
            tb[3,i,j] = data['b%dtoa%d' % (i+1,j+1)][0,3,0]

# write to file
    outfile = '.'.join(inpfile.split('.')[:-1]) + '.out'
    with open(outfile, 'w') as file:
        file.write('# Brightness temperatures of input model %s - model 0\n' % inpfile)
        file.write('%12s' % '# Freq (GHz)')
        for i in range(num_dirs):
            file.write('%10.1f' % amu[i])
        file.write('\n')
        for i in range(num_bands):
            file.write('%12.1f' % freq[i])
            for j in range(num_dirs):
                file.write('%10.2f' % tb[0,i,j])
            file.write('\n')

        file.write('# Brightness temperatures of input model %s - model 1\n' % inpfile)
        file.write('%12s' % '# Freq (GHz)')
        for i in range(num_dirs):
            file.write('%10.1f' % amu[i])
        file.write('\n')
        for i in range(num_bands):
            file.write('%12.1f' % freq[i])
            for j in range(num_dirs):
                file.write('%10.2f' % tb[1,i,j])
            file.write('\n')

        file.write('# Brightness temperatures of input model %s - model 2\n' % inpfile)
        file.write('%12s' % '# Freq (GHz)')
        for i in range(num_dirs):
            file.write('%10.1f' % amu[i])
        file.write('\n')
        for i in range(num_bands):
            file.write('%12.1f' % freq[i])
            for j in range(num_dirs):
                file.write('%10.2f' % tb[2,i,j])
            file.write('\n')

        file.write('# Brightness temperatures of input model %s - model 3\n' % inpfile)
        file.write('%12s' % '# Freq (GHz)')
        for i in range(num_dirs):
            file.write('%10.1f' % amu[i])
        file.write('\n')
        for i in range(num_bands):
            file.write('%12.1f' % freq[i])
            for j in range(num_dirs):
                file.write('%10.2f' % tb[3,i,j])
            file.write('\n')
    print('Brightness temperatures written to %s' % outfile)
    return outfile
