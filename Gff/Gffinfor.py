#! /usr/bin/env python3

__author__ = 'd2jvkpn'
__version__ = '2.1.2'
__release__ = '2018-07-24'
__project__ = 'https://github.com/d2jvkpn/dksh'
__license__ = 'GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)'

import argparse, os, getpass, time, string, subprocess
parser = argparse.ArgumentParser()

parser.add_argument('ImagePlus', nargs = '+', help='set local docker image ' + \
  ' to creat container, you have three work modes to choose, ' + \
  'mode1: by explicitly set user as empty (-u "") to create container and ' + \
    'enter shell environment (to explory and maintance image); ' + \
  'mode2: excute a command in conatainer ' + \
    '(mount host work directory in container); ' + \
  'mode3: run a script in container enviroment, a log file will be created ' + \
    'to save parameters, stdout, stderr, exit status, etc.')

parser.add_argument('-A', dest='A', default='', help='quoted string, ' + \
  'addition argument for docker run, except -w , --rm and --name.')

parser.add_argument('--mount', dest='mount', nargs='*', default=[], 
  help='path mount to container, e.g. /home/user/download /opt:/mnt/opt ../::wr' )

parser.add_argument('-w', dest='w', default=os.getcwd(),
  help='set work path of hostmachine, default: $PWD.')

parser.add_argument('-u', dest='u', default=getpass.getuser(),
  help='set login user in container, default: $USER.')

parser.add_argument('-n', dest='n', default=None, 
  help='set container name.')

parser.add_argument('-p', dest='p', type=int, default=None,
  help='maximum number of processes container can use.')

parser.add_argument('-m', dest='m', default=None, help='amount of ' +  \
  'memory(GB) container can use, followed by a suffix of b, k, m, g.')

parser.add_argument('-cpu_period', dest='cpu_period', type=int, default=100000, 
  help='set cpu_period, default: %(default)s.')

parser.add_argument('-shell', dest='shell', default='/bin/bash', 
  help='set container shell, default: %(default)s.')

if len(os.sys.argv) == 1 or os.sys.argv[1] in ['-h', '--help']:
    parser.print_help()

    print('\nNote: stop running script by kill container is recommanded.')

    _ = '\nauthor: {}\nversion: {}\nrelease: {}\nproject: {}\nlicense: {}\n'
    __ = [__author__,  __version__, __release__, __project__, __license__]
    print (_.format (*__))

    os.sys.exit(0)

args = parser.parse_args()

##
def validName (k):
    _ = list (string.ascii_letters + string.digits + '.-_')
    return(False not in [ i in _ for i in k ] )

def pathMount(path):
    sp = path.split(':')
    sp[0] = os.path.realpath(sp[0])

    _  = 'Error: mount path "{}" is not exist.'.format(sp[0])
    if not os.path.isdir(sp[0]): os.sys.exit(_)
    if path == '': print('Error: empty string for "--mount".')
    if len(sp) > 3: os.sys.exit('Error: wrong arugment "%s" for --mount' % path)

    if len(sp) == 2: sp.append('ro')
    if len(sp) == 1: sp += [os.path.abspath(sp[0]), 'ro']
    if sp[1] == '': sp[1] = os.path.abspath(sp[0])
    if sp[2] == '': sp[2] = 'ro'

    return('-v ' + ':'.join(sp))


####
_ = 'Error: only one or two position argument allowed.'
if len(args.ImagePlus) not in [1,2]: os.sys.exit(_)

image = args.ImagePlus[0]

####
err, output = subprocess.getstatusoutput('''docker images | wc -l''')

if err != 0: os.sys.exit('Error: failed to run docker command.')
if output == '1':
   print("Error: no docker image available.")
   os.sys.exit(1)

err, output = subprocess.getstatusoutput('''docker images |
awk -v OFS="\t" \'NR>1{if(NF==8) {s=$7" "$8} else {s=$8}; \
print $1":"$2, $3, s"\t"$4" "$5" "$6}\'''')

imageName = {}; imageID = {}; ilist = []

for i in output.split('\n'):
   ii = i.split('\t', 2)
   ilist.append(ii[0])
   imageName[ii[0]], imageID[ii[1]] = [ii[1], ii[2]], [ii[0], ii[2]]

if image in imageName:
    Image = image + ', ' + imageName[image][0]
elif image + ':latest' in imageName:
    Image = image + ':latest, ' + imageName[image + ':latest'][0]
elif image in imageID:
    Image = imageID[image][0] + ', ' + image
else:
    print('Cant\'t find "%s", available local image(s):' % image)
    _ = ['    ' + k + '\t' + '\t'.join(imageName[k]) for k in ilist]
    os.system("echo '%s' | column -t -s $'\t'" % '\n'.join(_) )
    os.sys.exit(1)


####
t0 = time.time()
t1 = time.strftime('%Y %m %d %H %M %S %s %z').split()
StartAt = '-'.join(t1[0:3]) + ' ' + ':'.join(t1[3:6]) + ' ' + t1[7]

if args.n == None:
    args.n =  hex(int(''.join(t1[:-1]).replace('.', ''))).replace('0x', '')

_ = 'Error: invalid container name "%s"'
if not validName (args.n): os.sys.exit(_ % args.n)

if args.u == "":
    print('Creating container "%s" from "%s"...' % (args.n, image))

    _ = 'docker run --name=%s -e PoweredBy=dksh --rm -it %s %s'
    os.system (_ % (args.n, image, args.shell))
    os.sys.exit(0)

##
mountPath = ' '.join ([pathMount(i) for i in args.mount])
HostID = subprocess.getoutput("hostid")
HostPID = os.getpid()

_  = 'Error: work path "%s" is not exist.' % args.w
if not os.path.isdir(args.w): os.sys.exit(_)

setWorkpath = '-v {0}:{1}'.format (os.path.realpath (args.w), '/mnt/HostPath')

_ = '--cpu-period %s --cpu-quota {}' % args.cpu_period
setCPU = '' if args.p == None else _ .format(args.cpu_period * args.p)

setMemory = '' if args.m == None else '-m %s' % args.m

##
cmd1 = ' '.join(['docker run %s --rm' % args.A, setWorkpath, 
    '-e DockerImage="%s"' % Image,
    '-e StartAt="%s"' % StartAt,
    '-e HostUser="%s" -e HostPath="%s"' % (args.u, os.path.abspath(args.w)), 
    '-e HostID="%s" -e HostPID=%s -e PoweredBy=dksh' % (HostID, HostPID),
    '-e Container=%s' % args.n, mountPath, setCPU, setMemory])

cmd2 = '{} {} -c "useradd {} &> /dev/null'.format (image, args.shell, args.u)

####
if len(args.ImagePlus) == 2:
    isScript = args.ImagePlus[1].endswith('.sh') and ' ' not in args.ImagePlus[1]
else:
    isScript = False

if not isScript:
    if len(args.ImagePlus) == 1:
        cmd2 += '; cd /mnt/HostPath; su %s -s %s"' % (args.u, args.shell)
    else:
        _ = '; su -l %s -c \'cd /mnt/HostPath && %s\'"'
        cmd2 += _ % (args.u, args.ImagePlus[1])

    _ = ' '.join ([cmd1, '--name=%s -it' % args.n, cmd2])

    print('Creating container "%s" from "%s"...' % (args.n, image))

    os.system (' '.join ([cmd1, '--name=%s -it' % args.n, cmd2]))
    os.sys.exit (0)

s = args.ImagePlus[1]

if not os.path.isfile (s): os.sys.exit ('Error: can\'t find script "%s".' % s)

if os.path.islink(s):
    cmd1 += ' ' + pathMount (os.path.dirname (os.path.realpath(s)))

script = os.path.abspath (os.path.dirname (s)) + '/' + os.path.basename(s)
cmd1 += ' ' + pathMount (os.path.dirname (script))

prefix = args.w + '/log/' + ''.join(t1[0:-2]) + '_' + args.n

os.system('mkdir -p %s/log' % args.w)

with open(prefix + '.logging', 'w') as f:
    f.write('\n## '.join(['## StartAt: %s' % StartAt,
    'HostUser: %s' % args.u,
    'HostPath: %s' % os.path.abspath (args.w),
    'HostID: %s' % HostID,
    'HostPID: %s' % HostPID,
    'Script: %s' % script,
    'DockerImage: %s' % Image,
    'Container: %s' % args.n,
    'PoweredBy: %s' % __project__,
    'CPU: %s' % ('' if args.p == None else args.p),
    'Memory: %s' % ('' if args.m == None else args.m),
    'DockerArgs: %s' % args.A,
    'Mount: %s\n' % '\n##   '.join(args.mount)]))

    f.write('\n\n')

_ = '; su -l %s -c \'cd /mnt/HostPath && %s %s\'" &>> %s.logging'
cmd2 += _ % (args.u, args.shell, script, prefix)
command = ' '.join([cmd1, '--name=%s' % args.n, cmd2])

_ = '\n## Running "%s" in "%s" ("%s")\n   WorkPath: %s' % \
(args.ImagePlus[1], args.n, image, args.w)

print (_ + '\n   SatrtAt: %s\n   PID: %s' % (StartAt, HostPID), \
file = open('run.log', 'a'))

try:    
    exitcode, _ = subprocess.getstatusoutput (command)
except KeyboardInterrupt:
    exitcode = 'interrupted'

_, __ = divmod (time.time () - t0, 3600)
elapsed = "%dd %dh %dm %.3fs" % (*divmod (_, 24), *divmod (__, 60))

with open (prefix + '.logging', 'a') as f:
    f.write ('\n## '.join (['\n', 
      'ExitCode: %s' % exitcode,
      'Elapsed: %s' % elapsed,
      'EndAt: %s\n' % time.strftime('%Y-%m-%d %H:%M:%S %z')])
    )

if isinstance (exitcode, int):
    suffix = 'log' if exitcode == 0 else 'failed'
else:
    suffix = 'interrupted'

os.system ('mv {0}.logging {0}.{1}'.format(prefix, suffix))
