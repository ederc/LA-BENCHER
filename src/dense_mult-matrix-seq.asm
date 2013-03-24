
dense_mult-matrix-seq.o:     file format elf64-x86-64


Disassembly of section .text:

0000000000000000 <_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.part.1>:
   0:	48 83 ec 08          	sub    $0x8,%rsp
   4:	48 8b 07             	mov    (%rdi),%rax
   7:	48 03 78 e8          	add    -0x18(%rax),%rdi
   b:	8b 77 20             	mov    0x20(%rdi),%esi
   e:	83 ce 01             	or     $0x1,%esi
  11:	e8 00 00 00 00       	callq  16 <_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.part.1+0x16>
  16:	48 83 c4 08          	add    $0x8,%rsp
  1a:	c3                   	retq   
  1b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000000020 <_Z7multSEQR6MatrixRKS_S2_ii>:
  20:	41 57                	push   %r15
  22:	41 56                	push   %r14
  24:	41 55                	push   %r13
  26:	41 54                	push   %r12
  28:	55                   	push   %rbp
  29:	53                   	push   %rbx
  2a:	48 81 ec e8 00 00 00 	sub    $0xe8,%rsp
  31:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
  38:	00 00 
  3a:	48 89 84 24 d8 00 00 	mov    %rax,0xd8(%rsp)
  41:	00 
  42:	31 c0                	xor    %eax,%eax
  44:	41 83 f8 01          	cmp    $0x1,%r8d
  48:	48 89 7c 24 40       	mov    %rdi,0x40(%rsp)
  4d:	48 89 74 24 58       	mov    %rsi,0x58(%rsp)
  52:	48 89 54 24 50       	mov    %rdx,0x50(%rsp)
  57:	89 4c 24 7c          	mov    %ecx,0x7c(%rsp)
  5b:	44 89 44 24 78       	mov    %r8d,0x78(%rsp)
  60:	0f 84 39 0f 00 00    	je     f9f <_Z7multSEQR6MatrixRKS_S2_ii+0xf7f>
  66:	48 8b 4c 24 58       	mov    0x58(%rsp),%rcx
  6b:	48 8b 54 24 50       	mov    0x50(%rsp),%rdx
  70:	8b 01                	mov    (%rcx),%eax
  72:	44 8b 7a 04          	mov    0x4(%rdx),%r15d
  76:	44 8b 2a             	mov    (%rdx),%r13d
  79:	89 44 24 60          	mov    %eax,0x60(%rsp)
  7d:	4c 8b 44 24 40       	mov    0x40(%rsp),%r8
  82:	8b 6c 24 60          	mov    0x60(%rsp),%ebp
  86:	4d 8b 48 10          	mov    0x10(%r8),%r9
  8a:	4d 8b 50 18          	mov    0x18(%r8),%r10
  8e:	41 0f af ef          	imul   %r15d,%ebp
  92:	4d 29 ca             	sub    %r9,%r10
  95:	49 c1 fa 02          	sar    $0x2,%r10
  99:	4c 39 d5             	cmp    %r10,%rbp
  9c:	0f 87 e6 0e 00 00    	ja     f88 <_Z7multSEQR6MatrixRKS_S2_ii+0xf68>
  a2:	73 0e                	jae    b2 <_Z7multSEQR6MatrixRKS_S2_ii+0x92>
  a4:	4c 8b 64 24 40       	mov    0x40(%rsp),%r12
  a9:	4d 8d 1c a9          	lea    (%r9,%rbp,4),%r11
  ad:	4d 89 5c 24 18       	mov    %r11,0x18(%r12)
  b2:	ba 15 00 00 00       	mov    $0x15,%edx
  b7:	be 00 00 00 00       	mov    $0x0,%esi
  bc:	bf 00 00 00 00       	mov    $0x0,%edi
  c1:	e8 00 00 00 00       	callq  c6 <_Z7multSEQR6MatrixRKS_S2_ii+0xa6>
  c6:	4c 8b 35 00 00 00 00 	mov    0x0(%rip),%r14        # cd <_Z7multSEQR6MatrixRKS_S2_ii+0xad>
  cd:	49 8b 4e e8          	mov    -0x18(%r14),%rcx
  d1:	48 8b 99 00 00 00 00 	mov    0x0(%rcx),%rbx
  d8:	48 85 db             	test   %rbx,%rbx
  db:	0f 84 f6 0e 00 00    	je     fd7 <_Z7multSEQR6MatrixRKS_S2_ii+0xfb7>
  e1:	80 7b 38 00          	cmpb   $0x0,0x38(%rbx)
  e5:	0f 84 f5 08 00 00    	je     9e0 <_Z7multSEQR6MatrixRKS_S2_ii+0x9c0>
  eb:	0f b6 43 43          	movzbl 0x43(%rbx),%eax
  ef:	0f be f0             	movsbl %al,%esi
  f2:	bf 00 00 00 00       	mov    $0x0,%edi
  f7:	e8 00 00 00 00       	callq  fc <_Z7multSEQR6MatrixRKS_S2_ii+0xdc>
  fc:	48 89 c7             	mov    %rax,%rdi
  ff:	e8 00 00 00 00       	callq  104 <_Z7multSEQR6MatrixRKS_S2_ii+0xe4>
 104:	48 8d bc 24 80 00 00 	lea    0x80(%rsp),%rdi
 10b:	00 
 10c:	31 f6                	xor    %esi,%esi
 10e:	e8 00 00 00 00       	callq  113 <_Z7multSEQR6MatrixRKS_S2_ii+0xf3>
 113:	e8 00 00 00 00       	callq  118 <_Z7multSEQR6MatrixRKS_S2_ii+0xf8>
 118:	83 7c 24 78 01       	cmpl   $0x1,0x78(%rsp)
 11d:	48 89 44 24 70       	mov    %rax,0x70(%rsp)
 122:	0f 84 b0 0a 00 00    	je     bd8 <_Z7multSEQR6MatrixRKS_S2_ii+0xbb8>
 128:	8b 74 24 60          	mov    0x60(%rsp),%esi
 12c:	45 31 db             	xor    %r11d,%r11d
 12f:	31 db                	xor    %ebx,%ebx
 131:	31 d2                	xor    %edx,%edx
 133:	0f 57 d2             	xorps  %xmm2,%xmm2
 136:	85 f6                	test   %esi,%esi
 138:	0f 84 87 02 00 00    	je     3c5 <_Z7multSEQR6MatrixRKS_S2_ii+0x3a5>
 13e:	4c 8b 54 24 58       	mov    0x58(%rsp),%r10
 143:	48 8b 6c 24 50       	mov    0x50(%rsp),%rbp
 148:	89 54 24 38          	mov    %edx,0x38(%rsp)
 14c:	0f 1f 40 00          	nopl   0x0(%rax)
 150:	45 85 ff             	test   %r15d,%r15d
 153:	47 8d 04 2b          	lea    (%r11,%r13,1),%r8d
 157:	0f 84 4d 02 00 00    	je     3aa <_Z7multSEQR6MatrixRKS_S2_ii+0x38a>
 15d:	48 8b 7c 24 40       	mov    0x40(%rsp),%rdi
 162:	47 8d 04 2b          	lea    (%r11,%r13,1),%r8d
 166:	4c 8b 4f 10          	mov    0x10(%rdi),%r9
 16a:	31 ff                	xor    %edi,%edi
 16c:	0f 1f 40 00          	nopl   0x0(%rax)
 170:	45 85 ed             	test   %r13d,%r13d
 173:	0f 84 5f 08 00 00    	je     9d8 <_Z7multSEQR6MatrixRKS_S2_ii+0x9b8>
 179:	49 8b 72 10          	mov    0x10(%r10),%rsi
 17d:	48 8b 4d 10          	mov    0x10(%rbp),%rcx
 181:	45 89 de             	mov    %r11d,%r14d
 184:	89 f8                	mov    %edi,%eax
 186:	45 89 dc             	mov    %r11d,%r12d
 189:	42 8d 14 3f          	lea    (%rdi,%r15,1),%edx
 18d:	41 f7 d4             	not    %r12d
 190:	f3 42 0f 10 04 b6    	movss  (%rsi,%r14,4),%xmm0
 196:	45 01 c4             	add    %r8d,%r12d
 199:	f3 0f 59 04 81       	mulss  (%rcx,%rax,4),%xmm0
 19e:	41 8d 43 01          	lea    0x1(%r11),%eax
 1a2:	41 83 e4 07          	and    $0x7,%r12d
 1a6:	44 39 c0             	cmp    %r8d,%eax
 1a9:	f3 0f 58 c2          	addss  %xmm2,%xmm0
 1ad:	0f 84 e2 01 00 00    	je     395 <_Z7multSEQR6MatrixRKS_S2_ii+0x375>
 1b3:	45 85 e4             	test   %r12d,%r12d
 1b6:	0f 84 f7 00 00 00    	je     2b3 <_Z7multSEQR6MatrixRKS_S2_ii+0x293>
 1bc:	41 83 fc 01          	cmp    $0x1,%r12d
 1c0:	0f 84 c7 00 00 00    	je     28d <_Z7multSEQR6MatrixRKS_S2_ii+0x26d>
 1c6:	41 83 fc 02          	cmp    $0x2,%r12d
 1ca:	0f 84 a1 00 00 00    	je     271 <_Z7multSEQR6MatrixRKS_S2_ii+0x251>
 1d0:	41 83 fc 03          	cmp    $0x3,%r12d
 1d4:	74 7f                	je     255 <_Z7multSEQR6MatrixRKS_S2_ii+0x235>
 1d6:	41 83 fc 04          	cmp    $0x4,%r12d
 1da:	74 5d                	je     239 <_Z7multSEQR6MatrixRKS_S2_ii+0x219>
 1dc:	41 83 fc 05          	cmp    $0x5,%r12d
 1e0:	74 3b                	je     21d <_Z7multSEQR6MatrixRKS_S2_ii+0x1fd>
 1e2:	41 83 fc 06          	cmp    $0x6,%r12d
 1e6:	74 19                	je     201 <_Z7multSEQR6MatrixRKS_S2_ii+0x1e1>
 1e8:	41 89 d4             	mov    %edx,%r12d
 1eb:	f3 0f 10 0c 86       	movss  (%rsi,%rax,4),%xmm1
 1f0:	f3 42 0f 59 0c a1    	mulss  (%rcx,%r12,4),%xmm1
 1f6:	41 8d 43 02          	lea    0x2(%r11),%eax
 1fa:	44 01 fa             	add    %r15d,%edx
 1fd:	f3 0f 58 c1          	addss  %xmm1,%xmm0
 201:	41 89 c6             	mov    %eax,%r14d
 204:	41 89 d4             	mov    %edx,%r12d
 207:	83 c0 01             	add    $0x1,%eax
 20a:	f3 42 0f 10 1c b6    	movss  (%rsi,%r14,4),%xmm3
 210:	44 01 fa             	add    %r15d,%edx
 213:	f3 42 0f 59 1c a1    	mulss  (%rcx,%r12,4),%xmm3
 219:	f3 0f 58 c3          	addss  %xmm3,%xmm0
 21d:	41 89 c6             	mov    %eax,%r14d
 220:	41 89 d4             	mov    %edx,%r12d
 223:	83 c0 01             	add    $0x1,%eax
 226:	f3 42 0f 10 24 b6    	movss  (%rsi,%r14,4),%xmm4
 22c:	44 01 fa             	add    %r15d,%edx
 22f:	f3 42 0f 59 24 a1    	mulss  (%rcx,%r12,4),%xmm4
 235:	f3 0f 58 c4          	addss  %xmm4,%xmm0
 239:	41 89 c6             	mov    %eax,%r14d
 23c:	41 89 d4             	mov    %edx,%r12d
 23f:	83 c0 01             	add    $0x1,%eax
 242:	f3 42 0f 10 2c b6    	movss  (%rsi,%r14,4),%xmm5
 248:	44 01 fa             	add    %r15d,%edx
 24b:	f3 42 0f 59 2c a1    	mulss  (%rcx,%r12,4),%xmm5
 251:	f3 0f 58 c5          	addss  %xmm5,%xmm0
 255:	41 89 c6             	mov    %eax,%r14d
 258:	41 89 d4             	mov    %edx,%r12d
 25b:	83 c0 01             	add    $0x1,%eax
 25e:	f3 42 0f 10 34 b6    	movss  (%rsi,%r14,4),%xmm6
 264:	44 01 fa             	add    %r15d,%edx
 267:	f3 42 0f 59 34 a1    	mulss  (%rcx,%r12,4),%xmm6
 26d:	f3 0f 58 c6          	addss  %xmm6,%xmm0
 271:	41 89 c6             	mov    %eax,%r14d
 274:	41 89 d4             	mov    %edx,%r12d
 277:	83 c0 01             	add    $0x1,%eax
 27a:	f3 42 0f 10 3c b6    	movss  (%rsi,%r14,4),%xmm7
 280:	44 01 fa             	add    %r15d,%edx
 283:	f3 42 0f 59 3c a1    	mulss  (%rcx,%r12,4),%xmm7
 289:	f3 0f 58 c7          	addss  %xmm7,%xmm0
 28d:	41 89 c6             	mov    %eax,%r14d
 290:	41 89 d4             	mov    %edx,%r12d
 293:	83 c0 01             	add    $0x1,%eax
 296:	f3 46 0f 10 04 b6    	movss  (%rsi,%r14,4),%xmm8
 29c:	44 01 fa             	add    %r15d,%edx
 29f:	f3 46 0f 59 04 a1    	mulss  (%rcx,%r12,4),%xmm8
 2a5:	44 39 c0             	cmp    %r8d,%eax
 2a8:	f3 41 0f 58 c0       	addss  %xmm8,%xmm0
 2ad:	0f 84 e2 00 00 00    	je     395 <_Z7multSEQR6MatrixRKS_S2_ii+0x375>
 2b3:	41 89 c6             	mov    %eax,%r14d
 2b6:	41 89 d4             	mov    %edx,%r12d
 2b9:	44 01 fa             	add    %r15d,%edx
 2bc:	f3 46 0f 10 0c b6    	movss  (%rsi,%r14,4),%xmm9
 2c2:	44 8d 70 01          	lea    0x1(%rax),%r14d
 2c6:	f3 46 0f 59 0c a1    	mulss  (%rcx,%r12,4),%xmm9
 2cc:	41 89 d4             	mov    %edx,%r12d
 2cf:	44 01 fa             	add    %r15d,%edx
 2d2:	f3 46 0f 10 14 b6    	movss  (%rsi,%r14,4),%xmm10
 2d8:	44 8d 70 02          	lea    0x2(%rax),%r14d
 2dc:	f3 46 0f 59 14 a1    	mulss  (%rcx,%r12,4),%xmm10
 2e2:	41 89 d4             	mov    %edx,%r12d
 2e5:	44 01 fa             	add    %r15d,%edx
 2e8:	f3 46 0f 10 1c b6    	movss  (%rsi,%r14,4),%xmm11
 2ee:	44 8d 70 03          	lea    0x3(%rax),%r14d
 2f2:	f3 46 0f 59 1c a1    	mulss  (%rcx,%r12,4),%xmm11
 2f8:	41 89 d4             	mov    %edx,%r12d
 2fb:	44 01 fa             	add    %r15d,%edx
 2fe:	f3 46 0f 10 24 b6    	movss  (%rsi,%r14,4),%xmm12
 304:	44 8d 70 04          	lea    0x4(%rax),%r14d
 308:	f3 41 0f 58 c1       	addss  %xmm9,%xmm0
 30d:	f3 46 0f 59 24 a1    	mulss  (%rcx,%r12,4),%xmm12
 313:	41 89 d4             	mov    %edx,%r12d
 316:	f3 46 0f 10 2c b6    	movss  (%rsi,%r14,4),%xmm13
 31c:	44 8d 70 05          	lea    0x5(%rax),%r14d
 320:	f3 46 0f 59 2c a1    	mulss  (%rcx,%r12,4),%xmm13
 326:	44 01 fa             	add    %r15d,%edx
 329:	41 89 d4             	mov    %edx,%r12d
 32c:	f3 46 0f 10 34 b6    	movss  (%rsi,%r14,4),%xmm14
 332:	f3 41 0f 58 c2       	addss  %xmm10,%xmm0
 337:	f3 46 0f 59 34 a1    	mulss  (%rcx,%r12,4),%xmm14
 33d:	44 8d 70 06          	lea    0x6(%rax),%r14d
 341:	44 01 fa             	add    %r15d,%edx
 344:	41 89 d4             	mov    %edx,%r12d
 347:	f3 46 0f 10 3c b6    	movss  (%rsi,%r14,4),%xmm15
 34d:	f3 46 0f 59 3c a1    	mulss  (%rcx,%r12,4),%xmm15
 353:	44 8d 70 07          	lea    0x7(%rax),%r14d
 357:	44 01 fa             	add    %r15d,%edx
 35a:	f3 41 0f 58 c3       	addss  %xmm11,%xmm0
 35f:	41 89 d4             	mov    %edx,%r12d
 362:	83 c0 08             	add    $0x8,%eax
 365:	f3 42 0f 10 0c b6    	movss  (%rsi,%r14,4),%xmm1
 36b:	44 01 fa             	add    %r15d,%edx
 36e:	f3 42 0f 59 0c a1    	mulss  (%rcx,%r12,4),%xmm1
 374:	44 39 c0             	cmp    %r8d,%eax
 377:	f3 41 0f 58 c4       	addss  %xmm12,%xmm0
 37c:	f3 41 0f 58 c5       	addss  %xmm13,%xmm0
 381:	f3 41 0f 58 c6       	addss  %xmm14,%xmm0
 386:	f3 41 0f 58 c7       	addss  %xmm15,%xmm0
 38b:	f3 0f 58 c1          	addss  %xmm1,%xmm0
 38f:	0f 85 1e ff ff ff    	jne    2b3 <_Z7multSEQR6MatrixRKS_S2_ii+0x293>
 395:	8d 34 3b             	lea    (%rbx,%rdi,1),%esi
 398:	83 c7 01             	add    $0x1,%edi
 39b:	44 39 ff             	cmp    %r15d,%edi
 39e:	f3 41 0f 11 04 b1    	movss  %xmm0,(%r9,%rsi,4)
 3a4:	0f 85 c6 fd ff ff    	jne    170 <_Z7multSEQR6MatrixRKS_S2_ii+0x150>
 3aa:	83 44 24 38 01       	addl   $0x1,0x38(%rsp)
 3af:	44 01 fb             	add    %r15d,%ebx
 3b2:	44 8b 4c 24 60       	mov    0x60(%rsp),%r9d
 3b7:	44 39 4c 24 38       	cmp    %r9d,0x38(%rsp)
 3bc:	45 89 c3             	mov    %r8d,%r11d
 3bf:	0f 85 8b fd ff ff    	jne    150 <_Z7multSEQR6MatrixRKS_S2_ii+0x130>
 3c5:	48 8d bc 24 90 00 00 	lea    0x90(%rsp),%rdi
 3cc:	00 
 3cd:	31 f6                	xor    %esi,%esi
 3cf:	e8 00 00 00 00       	callq  3d4 <_Z7multSEQR6MatrixRKS_S2_ii+0x3b4>
 3d4:	e8 00 00 00 00       	callq  3d9 <_Z7multSEQR6MatrixRKS_S2_ii+0x3b9>
 3d9:	ba 33 00 00 00       	mov    $0x33,%edx
 3de:	be 00 00 00 00       	mov    $0x0,%esi
 3e3:	bf 00 00 00 00       	mov    $0x0,%edi
 3e8:	48 89 c3             	mov    %rax,%rbx
 3eb:	e8 00 00 00 00       	callq  3f0 <_Z7multSEQR6MatrixRKS_S2_ii+0x3d0>
 3f0:	48 8b 2d 00 00 00 00 	mov    0x0(%rip),%rbp        # 3f7 <_Z7multSEQR6MatrixRKS_S2_ii+0x3d7>
 3f7:	4c 8b 65 e8          	mov    -0x18(%rbp),%r12
 3fb:	4d 8b b4 24 00 00 00 	mov    0x0(%r12),%r14
 402:	00 
 403:	4d 85 f6             	test   %r14,%r14
 406:	0f 84 cb 0b 00 00    	je     fd7 <_Z7multSEQR6MatrixRKS_S2_ii+0xfb7>
 40c:	41 80 7e 38 00       	cmpb   $0x0,0x38(%r14)
 411:	0f 84 87 07 00 00    	je     b9e <_Z7multSEQR6MatrixRKS_S2_ii+0xb7e>
 417:	41 0f b6 46 43       	movzbl 0x43(%r14),%eax
 41c:	0f be f0             	movsbl %al,%esi
 41f:	bf 00 00 00 00       	mov    $0x0,%edi
 424:	e8 00 00 00 00       	callq  429 <_Z7multSEQR6MatrixRKS_S2_ii+0x409>
 429:	48 89 c7             	mov    %rax,%rdi
 42c:	e8 00 00 00 00       	callq  431 <_Z7multSEQR6MatrixRKS_S2_ii+0x411>
 431:	ba 20 00 00 00       	mov    $0x20,%edx
 436:	bf 00 00 00 00       	mov    $0x0,%edi
 43b:	be 00 00 00 00       	mov    $0x0,%esi
 440:	e8 00 00 00 00       	callq  445 <_Z7multSEQR6MatrixRKS_S2_ii+0x425>
 445:	48 8b 3d 00 00 00 00 	mov    0x0(%rip),%rdi        # 44c <_Z7multSEQR6MatrixRKS_S2_ii+0x42c>
 44c:	48 8b 57 e8          	mov    -0x18(%rdi),%rdx
 450:	4c 8b ba 00 00 00 00 	mov    0x0(%rdx),%r15
 457:	4d 85 ff             	test   %r15,%r15
 45a:	0f 84 77 0b 00 00    	je     fd7 <_Z7multSEQR6MatrixRKS_S2_ii+0xfb7>
 460:	41 80 7f 38 00       	cmpb   $0x0,0x38(%r15)
 465:	0f 84 18 07 00 00    	je     b83 <_Z7multSEQR6MatrixRKS_S2_ii+0xb63>
 46b:	41 0f b6 47 43       	movzbl 0x43(%r15),%eax
 470:	0f be f0             	movsbl %al,%esi
 473:	bf 00 00 00 00       	mov    $0x0,%edi
 478:	e8 00 00 00 00       	callq  47d <_Z7multSEQR6MatrixRKS_S2_ii+0x45d>
 47d:	48 89 c7             	mov    %rax,%rdi
 480:	e8 00 00 00 00       	callq  485 <_Z7multSEQR6MatrixRKS_S2_ii+0x465>
 485:	ba 12 00 00 00       	mov    $0x12,%edx
 48a:	be 00 00 00 00       	mov    $0x0,%esi
 48f:	bf 00 00 00 00       	mov    $0x0,%edi
 494:	e8 00 00 00 00       	callq  499 <_Z7multSEQR6MatrixRKS_S2_ii+0x479>
 499:	83 7c 24 78 01       	cmpl   $0x1,0x78(%rsp)
 49e:	ba 01 00 00 00       	mov    $0x1,%edx
 4a3:	be 00 00 00 00       	mov    $0x0,%esi
 4a8:	74 05                	je     4af <_Z7multSEQR6MatrixRKS_S2_ii+0x48f>
 4aa:	be 00 00 00 00       	mov    $0x0,%esi
 4af:	bf 00 00 00 00       	mov    $0x0,%edi
 4b4:	e8 00 00 00 00       	callq  4b9 <_Z7multSEQR6MatrixRKS_S2_ii+0x499>
 4b9:	4c 8b 1d 00 00 00 00 	mov    0x0(%rip),%r11        # 4c0 <_Z7multSEQR6MatrixRKS_S2_ii+0x4a0>
 4c0:	49 8b 73 e8          	mov    -0x18(%r11),%rsi
 4c4:	48 8b ae 00 00 00 00 	mov    0x0(%rsi),%rbp
 4cb:	48 85 ed             	test   %rbp,%rbp
 4ce:	0f 84 03 0b 00 00    	je     fd7 <_Z7multSEQR6MatrixRKS_S2_ii+0xfb7>
 4d4:	80 7d 38 00          	cmpb   $0x0,0x38(%rbp)
 4d8:	0f 84 dc 06 00 00    	je     bba <_Z7multSEQR6MatrixRKS_S2_ii+0xb9a>
 4de:	0f b6 45 43          	movzbl 0x43(%rbp),%eax
 4e2:	0f be f0             	movsbl %al,%esi
 4e5:	bf 00 00 00 00       	mov    $0x0,%edi
 4ea:	e8 00 00 00 00       	callq  4ef <_Z7multSEQR6MatrixRKS_S2_ii+0x4cf>
 4ef:	48 89 c7             	mov    %rax,%rdi
 4f2:	e8 00 00 00 00       	callq  4f7 <_Z7multSEQR6MatrixRKS_S2_ii+0x4d7>
 4f7:	48 8b 4c 24 58       	mov    0x58(%rsp),%rcx
 4fc:	4c 8b 54 24 50       	mov    0x50(%rsp),%r10
 501:	44 8b 09             	mov    (%rcx),%r9d
 504:	45 8b 02             	mov    (%r10),%r8d
 507:	45 8b 72 04          	mov    0x4(%r10),%r14d
 50b:	4d 01 c9             	add    %r9,%r9
 50e:	4d 0f af c8          	imul   %r8,%r9
 512:	4d 0f af ce          	imul   %r14,%r9
 516:	4d 85 c9             	test   %r9,%r9
 519:	0f 88 52 0a 00 00    	js     f71 <_Z7multSEQR6MatrixRKS_S2_ii+0xf51>
 51f:	f2 49 0f 2a e9       	cvtsi2sd %r9,%xmm5
 524:	f2 0f 11 6c 24 48    	movsd  %xmm5,0x48(%rsp)
 52a:	4c 8b ac 24 90 00 00 	mov    0x90(%rsp),%r13
 531:	00 
 532:	4c 2b ac 24 80 00 00 	sub    0x80(%rsp),%r13
 539:	00 
 53a:	b9 00 00 00 00       	mov    $0x0,%ecx
 53f:	48 8b bc 24 98 00 00 	mov    0x98(%rsp),%rdi
 546:	00 
 547:	48 2b bc 24 88 00 00 	sub    0x88(%rsp),%rdi
 54e:	00 
 54f:	ba 32 00 00 00       	mov    $0x32,%edx
 554:	48 2b 5c 24 70       	sub    0x70(%rsp),%rbx
 559:	f2 44 0f 10 05 00 00 	movsd  0x0(%rip),%xmm8        # 562 <_Z7multSEQR6MatrixRKS_S2_ii+0x542>
 560:	00 00 
 562:	be 01 00 00 00       	mov    $0x1,%esi
 567:	b8 01 00 00 00       	mov    $0x1,%eax
 56c:	f2 49 0f 2a fd       	cvtsi2sd %r13,%xmm7
 571:	f2 4c 0f 2a cf       	cvtsi2sd %rdi,%xmm9
 576:	48 8d bc 24 a0 00 00 	lea    0xa0(%rsp),%rdi
 57d:	00 
 57e:	f2 4c 0f 2a d3       	cvtsi2sd %rbx,%xmm10
 583:	f2 41 0f 59 f8       	mulsd  %xmm8,%xmm7
 588:	f2 45 0f 5e d0       	divsd  %xmm8,%xmm10
 58d:	f2 0f 11 7c 24 38    	movsd  %xmm7,0x38(%rsp)
 593:	f2 44 0f 58 4c 24 38 	addsd  0x38(%rsp),%xmm9
 59a:	f2 45 0f 5e c8       	divsd  %xmm8,%xmm9
 59f:	66 41 0f 28 c2       	movapd %xmm10,%xmm0
 5a4:	f2 44 0f 11 54 24 40 	movsd  %xmm10,0x40(%rsp)
 5ab:	f2 44 0f 11 4c 24 38 	movsd  %xmm9,0x38(%rsp)
 5b2:	e8 00 00 00 00       	callq  5b7 <_Z7multSEQR6MatrixRKS_S2_ii+0x597>
 5b7:	ba 12 00 00 00       	mov    $0x12,%edx
 5bc:	be 00 00 00 00       	mov    $0x0,%esi
 5c1:	bf 00 00 00 00       	mov    $0x0,%edi
 5c6:	41 89 c7             	mov    %eax,%r15d
 5c9:	e8 00 00 00 00       	callq  5ce <_Z7multSEQR6MatrixRKS_S2_ii+0x5ae>
 5ce:	be 01 00 00 00       	mov    $0x1,%esi
 5d3:	bf 00 00 00 00       	mov    $0x0,%edi
 5d8:	e8 00 00 00 00       	callq  5dd <_Z7multSEQR6MatrixRKS_S2_ii+0x5bd>
 5dd:	48 8b 18             	mov    (%rax),%rbx
 5e0:	48 89 c5             	mov    %rax,%rbp
 5e3:	48 8b 53 e8          	mov    -0x18(%rbx),%rdx
 5e7:	4c 8b a4 10 f0 00 00 	mov    0xf0(%rax,%rdx,1),%r12
 5ee:	00 
 5ef:	4d 85 e4             	test   %r12,%r12
 5f2:	0f 84 df 09 00 00    	je     fd7 <_Z7multSEQR6MatrixRKS_S2_ii+0xfb7>
 5f8:	41 80 7c 24 38 00    	cmpb   $0x0,0x38(%r12)
 5fe:	0f 84 63 05 00 00    	je     b67 <_Z7multSEQR6MatrixRKS_S2_ii+0xb47>
 604:	41 0f b6 44 24 43    	movzbl 0x43(%r12),%eax
 60a:	0f be f0             	movsbl %al,%esi
 60d:	48 89 ef             	mov    %rbp,%rdi
 610:	e8 00 00 00 00       	callq  615 <_Z7multSEQR6MatrixRKS_S2_ii+0x5f5>
 615:	48 89 c7             	mov    %rax,%rdi
 618:	e8 00 00 00 00       	callq  61d <_Z7multSEQR6MatrixRKS_S2_ii+0x5fd>
 61d:	ba 12 00 00 00       	mov    $0x12,%edx
 622:	be 00 00 00 00       	mov    $0x0,%esi
 627:	bf 00 00 00 00       	mov    $0x0,%edi
 62c:	e8 00 00 00 00       	callq  631 <_Z7multSEQR6MatrixRKS_S2_ii+0x611>
 631:	8b 74 24 7c          	mov    0x7c(%rsp),%esi
 635:	bf 00 00 00 00       	mov    $0x0,%edi
 63a:	e8 00 00 00 00       	callq  63f <_Z7multSEQR6MatrixRKS_S2_ii+0x61f>
 63f:	4c 8b 18             	mov    (%rax),%r11
 642:	49 89 c6             	mov    %rax,%r14
 645:	49 8b 73 e8          	mov    -0x18(%r11),%rsi
 649:	4c 8b ac 30 f0 00 00 	mov    0xf0(%rax,%rsi,1),%r13
 650:	00 
 651:	4d 85 ed             	test   %r13,%r13
 654:	0f 84 7d 09 00 00    	je     fd7 <_Z7multSEQR6MatrixRKS_S2_ii+0xfb7>
 65a:	41 80 7d 38 00       	cmpb   $0x0,0x38(%r13)
 65f:	0f 84 e6 04 00 00    	je     b4b <_Z7multSEQR6MatrixRKS_S2_ii+0xb2b>
 665:	41 0f b6 45 43       	movzbl 0x43(%r13),%eax
 66a:	0f be f0             	movsbl %al,%esi
 66d:	4c 89 f7             	mov    %r14,%rdi
 670:	e8 00 00 00 00       	callq  675 <_Z7multSEQR6MatrixRKS_S2_ii+0x655>
 675:	48 89 c7             	mov    %rax,%rdi
 678:	e8 00 00 00 00       	callq  67d <_Z7multSEQR6MatrixRKS_S2_ii+0x65d>
 67d:	ba 33 00 00 00       	mov    $0x33,%edx
 682:	be 00 00 00 00       	mov    $0x0,%esi
 687:	bf 00 00 00 00       	mov    $0x0,%edi
 68c:	e8 00 00 00 00       	callq  691 <_Z7multSEQR6MatrixRKS_S2_ii+0x671>
 691:	48 8b 0d 00 00 00 00 	mov    0x0(%rip),%rcx        # 698 <_Z7multSEQR6MatrixRKS_S2_ii+0x678>
 698:	4c 8b 49 e8          	mov    -0x18(%rcx),%r9
 69c:	49 8b 99 00 00 00 00 	mov    0x0(%r9),%rbx
 6a3:	48 85 db             	test   %rbx,%rbx
 6a6:	0f 84 2b 09 00 00    	je     fd7 <_Z7multSEQR6MatrixRKS_S2_ii+0xfb7>
 6ac:	80 7b 38 00          	cmpb   $0x0,0x38(%rbx)
 6b0:	0f 84 78 04 00 00    	je     b2e <_Z7multSEQR6MatrixRKS_S2_ii+0xb0e>
 6b6:	0f b6 43 43          	movzbl 0x43(%rbx),%eax
 6ba:	0f be f0             	movsbl %al,%esi
 6bd:	bf 00 00 00 00       	mov    $0x0,%edi
 6c2:	e8 00 00 00 00       	callq  6c7 <_Z7multSEQR6MatrixRKS_S2_ii+0x6a7>
 6c7:	48 89 c7             	mov    %rax,%rdi
 6ca:	e8 00 00 00 00       	callq  6cf <_Z7multSEQR6MatrixRKS_S2_ii+0x6af>
 6cf:	be 00 00 00 00       	mov    $0x0,%esi
 6d4:	ba 12 00 00 00       	mov    $0x12,%edx
 6d9:	bf 00 00 00 00       	mov    $0x0,%edi
 6de:	e8 00 00 00 00       	callq  6e3 <_Z7multSEQR6MatrixRKS_S2_ii+0x6c3>
 6e3:	4c 8b 05 00 00 00 00 	mov    0x0(%rip),%r8        # 6ea <_Z7multSEQR6MatrixRKS_S2_ii+0x6ca>
 6ea:	45 8d 57 05          	lea    0x5(%r15),%r10d
 6ee:	f2 0f 10 44 24 38    	movsd  0x38(%rsp),%xmm0
 6f4:	4d 63 fa             	movslq %r10d,%r15
 6f7:	49 8b 78 e8          	mov    -0x18(%r8),%rdi
 6fb:	4c 89 bf 00 00 00 00 	mov    %r15,0x0(%rdi)
 702:	49 8b 50 e8          	mov    -0x18(%r8),%rdx
 706:	bf 00 00 00 00       	mov    $0x0,%edi
 70b:	48 c7 82 00 00 00 00 	movq   $0x4,0x0(%rdx)
 712:	04 00 00 00 
 716:	4d 8b 70 e8          	mov    -0x18(%r8),%r14
 71a:	41 8b 86 00 00 00 00 	mov    0x0(%r14),%eax
 721:	25 fb fe ff ff       	and    $0xfffffefb,%eax
 726:	83 c8 04             	or     $0x4,%eax
 729:	41 89 86 00 00 00 00 	mov    %eax,0x0(%r14)
 730:	e8 00 00 00 00       	callq  735 <_Z7multSEQR6MatrixRKS_S2_ii+0x715>
 735:	be 00 00 00 00       	mov    $0x0,%esi
 73a:	49 89 c5             	mov    %rax,%r13
 73d:	ba 04 00 00 00       	mov    $0x4,%edx
 742:	48 89 c7             	mov    %rax,%rdi
 745:	e8 00 00 00 00       	callq  74a <_Z7multSEQR6MatrixRKS_S2_ii+0x72a>
 74a:	4d 8b 5d 00          	mov    0x0(%r13),%r11
 74e:	49 8b 73 e8          	mov    -0x18(%r11),%rsi
 752:	49 8b ac 35 f0 00 00 	mov    0xf0(%r13,%rsi,1),%rbp
 759:	00 
 75a:	48 85 ed             	test   %rbp,%rbp
 75d:	0f 84 74 08 00 00    	je     fd7 <_Z7multSEQR6MatrixRKS_S2_ii+0xfb7>
 763:	80 7d 38 00          	cmpb   $0x0,0x38(%rbp)
 767:	0f 84 a5 03 00 00    	je     b12 <_Z7multSEQR6MatrixRKS_S2_ii+0xaf2>
 76d:	0f b6 45 43          	movzbl 0x43(%rbp),%eax
 771:	0f be f0             	movsbl %al,%esi
 774:	4c 89 ef             	mov    %r13,%rdi
 777:	e8 00 00 00 00       	callq  77c <_Z7multSEQR6MatrixRKS_S2_ii+0x75c>
 77c:	48 89 c7             	mov    %rax,%rdi
 77f:	e8 00 00 00 00       	callq  784 <_Z7multSEQR6MatrixRKS_S2_ii+0x764>
 784:	ba 12 00 00 00       	mov    $0x12,%edx
 789:	be 00 00 00 00       	mov    $0x0,%esi
 78e:	bf 00 00 00 00       	mov    $0x0,%edi
 793:	e8 00 00 00 00       	callq  798 <_Z7multSEQR6MatrixRKS_S2_ii+0x778>
 798:	48 8b 0d 00 00 00 00 	mov    0x0(%rip),%rcx        # 79f <_Z7multSEQR6MatrixRKS_S2_ii+0x77f>
 79f:	f2 0f 10 44 24 40    	movsd  0x40(%rsp),%xmm0
 7a5:	bf 00 00 00 00       	mov    $0x0,%edi
 7aa:	4c 8b 49 e8          	mov    -0x18(%rcx),%r9
 7ae:	4d 89 b9 00 00 00 00 	mov    %r15,0x0(%r9)
 7b5:	4c 8b 61 e8          	mov    -0x18(%rcx),%r12
 7b9:	49 c7 84 24 00 00 00 	movq   $0x4,0x0(%r12)
 7c0:	00 04 00 00 00 
 7c5:	4c 8b 51 e8          	mov    -0x18(%rcx),%r10
 7c9:	45 8b 82 00 00 00 00 	mov    0x0(%r10),%r8d
 7d0:	41 81 e0 fb fe ff ff 	and    $0xfffffefb,%r8d
 7d7:	41 83 c8 04          	or     $0x4,%r8d
 7db:	45 89 82 00 00 00 00 	mov    %r8d,0x0(%r10)
 7e2:	e8 00 00 00 00       	callq  7e7 <_Z7multSEQR6MatrixRKS_S2_ii+0x7c7>
 7e7:	ba 04 00 00 00       	mov    $0x4,%edx
 7ec:	49 89 c6             	mov    %rax,%r14
 7ef:	48 89 c7             	mov    %rax,%rdi
 7f2:	be 00 00 00 00       	mov    $0x0,%esi
 7f7:	e8 00 00 00 00       	callq  7fc <_Z7multSEQR6MatrixRKS_S2_ii+0x7dc>
 7fc:	49 8b 3e             	mov    (%r14),%rdi
 7ff:	48 8b 57 e8          	mov    -0x18(%rdi),%rdx
 803:	4d 8b ac 16 f0 00 00 	mov    0xf0(%r14,%rdx,1),%r13
 80a:	00 
 80b:	4d 85 ed             	test   %r13,%r13
 80e:	0f 84 c3 07 00 00    	je     fd7 <_Z7multSEQR6MatrixRKS_S2_ii+0xfb7>
 814:	41 80 7d 38 00       	cmpb   $0x0,0x38(%r13)
 819:	0f 84 d7 02 00 00    	je     af6 <_Z7multSEQR6MatrixRKS_S2_ii+0xad6>
 81f:	41 0f b6 45 43       	movzbl 0x43(%r13),%eax
 824:	0f be f0             	movsbl %al,%esi
 827:	4c 89 f7             	mov    %r14,%rdi
 82a:	e8 00 00 00 00       	callq  82f <_Z7multSEQR6MatrixRKS_S2_ii+0x80f>
 82f:	48 89 c7             	mov    %rax,%rdi
 832:	e8 00 00 00 00       	callq  837 <_Z7multSEQR6MatrixRKS_S2_ii+0x817>
 837:	f2 44 0f 10 5c 24 40 	movsd  0x40(%rsp),%xmm11
 83e:	66 44 0f 2e 1d 00 00 	ucomisd 0x0(%rip),%xmm11        # 847 <_Z7multSEQR6MatrixRKS_S2_ii+0x827>
 845:	00 00 
 847:	0f 87 02 02 00 00    	ja     a4f <_Z7multSEQR6MatrixRKS_S2_ii+0xa2f>
 84d:	bf 00 00 00 00       	mov    $0x0,%edi
 852:	ba 33 00 00 00       	mov    $0x33,%edx
 857:	be 00 00 00 00       	mov    $0x0,%esi
 85c:	e8 00 00 00 00       	callq  861 <_Z7multSEQR6MatrixRKS_S2_ii+0x841>
 861:	4c 8b 05 00 00 00 00 	mov    0x0(%rip),%r8        # 868 <_Z7multSEQR6MatrixRKS_S2_ii+0x848>
 868:	49 8b 78 e8          	mov    -0x18(%r8),%rdi
 86c:	48 8b af 00 00 00 00 	mov    0x0(%rdi),%rbp
 873:	48 85 ed             	test   %rbp,%rbp
 876:	0f 84 5b 07 00 00    	je     fd7 <_Z7multSEQR6MatrixRKS_S2_ii+0xfb7>
 87c:	80 7d 38 00          	cmpb   $0x0,0x38(%rbp)
 880:	0f 84 ad 01 00 00    	je     a33 <_Z7multSEQR6MatrixRKS_S2_ii+0xa13>
 886:	0f b6 45 43          	movzbl 0x43(%rbp),%eax
 88a:	0f be f0             	movsbl %al,%esi
 88d:	bf 00 00 00 00       	mov    $0x0,%edi
 892:	e8 00 00 00 00       	callq  897 <_Z7multSEQR6MatrixRKS_S2_ii+0x877>
 897:	48 89 c7             	mov    %rax,%rdi
 89a:	e8 00 00 00 00       	callq  89f <_Z7multSEQR6MatrixRKS_S2_ii+0x87f>
 89f:	f2 44 0f 10 64 24 38 	movsd  0x38(%rsp),%xmm12
 8a6:	ba 12 00 00 00       	mov    $0x12,%edx
 8ab:	f2 44 0f 10 6c 24 48 	movsd  0x48(%rsp),%xmm13
 8b2:	be 00 00 00 00       	mov    $0x0,%esi
 8b7:	f2 44 0f 59 25 00 00 	mulsd  0x0(%rip),%xmm12        # 8c0 <_Z7multSEQR6MatrixRKS_S2_ii+0x8a0>
 8be:	00 00 
 8c0:	bf 00 00 00 00       	mov    $0x0,%edi
 8c5:	f2 45 0f 5e ec       	divsd  %xmm12,%xmm13
 8ca:	f2 44 0f 11 6c 24 20 	movsd  %xmm13,0x20(%rsp)
 8d1:	e8 00 00 00 00       	callq  8d6 <_Z7multSEQR6MatrixRKS_S2_ii+0x8b6>
 8d6:	4c 8b 1d 00 00 00 00 	mov    0x0(%rip),%r11        # 8dd <_Z7multSEQR6MatrixRKS_S2_ii+0x8bd>
 8dd:	bf 00 00 00 00       	mov    $0x0,%edi
 8e2:	f2 0f 10 44 24 20    	movsd  0x20(%rsp),%xmm0
 8e8:	49 8b 53 e8          	mov    -0x18(%r11),%rdx
 8ec:	4c 89 ba 00 00 00 00 	mov    %r15,0x0(%rdx)
 8f3:	4d 8b 7b e8          	mov    -0x18(%r11),%r15
 8f7:	49 c7 87 00 00 00 00 	movq   $0x4,0x0(%r15)
 8fe:	04 00 00 00 
 902:	49 8b 73 e8          	mov    -0x18(%r11),%rsi
 906:	8b 9e 00 00 00 00    	mov    0x0(%rsi),%ebx
 90c:	81 e3 fb fe ff ff    	and    $0xfffffefb,%ebx
 912:	83 cb 04             	or     $0x4,%ebx
 915:	89 9e 00 00 00 00    	mov    %ebx,0x0(%rsi)
 91b:	e8 00 00 00 00       	callq  920 <_Z7multSEQR6MatrixRKS_S2_ii+0x900>
 920:	48 8b 08             	mov    (%rax),%rcx
 923:	49 89 c4             	mov    %rax,%r12
 926:	4c 8b 49 e8          	mov    -0x18(%rcx),%r9
 92a:	4e 8b b4 08 f0 00 00 	mov    0xf0(%rax,%r9,1),%r14
 931:	00 
 932:	4d 85 f6             	test   %r14,%r14
 935:	0f 84 9c 06 00 00    	je     fd7 <_Z7multSEQR6MatrixRKS_S2_ii+0xfb7>
 93b:	41 80 7e 38 00       	cmpb   $0x0,0x38(%r14)
 940:	0f 84 d1 00 00 00    	je     a17 <_Z7multSEQR6MatrixRKS_S2_ii+0x9f7>
 946:	41 0f b6 46 43       	movzbl 0x43(%r14),%eax
 94b:	0f be f0             	movsbl %al,%esi
 94e:	4c 89 e7             	mov    %r12,%rdi
 951:	e8 00 00 00 00       	callq  956 <_Z7multSEQR6MatrixRKS_S2_ii+0x936>
 956:	48 89 c7             	mov    %rax,%rdi
 959:	e8 00 00 00 00       	callq  95e <_Z7multSEQR6MatrixRKS_S2_ii+0x93e>
 95e:	ba 33 00 00 00       	mov    $0x33,%edx
 963:	be 00 00 00 00       	mov    $0x0,%esi
 968:	bf 00 00 00 00       	mov    $0x0,%edi
 96d:	e8 00 00 00 00       	callq  972 <_Z7multSEQR6MatrixRKS_S2_ii+0x952>
 972:	4c 8b 15 00 00 00 00 	mov    0x0(%rip),%r10        # 979 <_Z7multSEQR6MatrixRKS_S2_ii+0x959>
 979:	4d 8b 42 e8          	mov    -0x18(%r10),%r8
 97d:	49 8b a8 00 00 00 00 	mov    0x0(%r8),%rbp
 984:	48 85 ed             	test   %rbp,%rbp
 987:	0f 84 4a 06 00 00    	je     fd7 <_Z7multSEQR6MatrixRKS_S2_ii+0xfb7>
 98d:	80 7d 38 00          	cmpb   $0x0,0x38(%rbp)
 991:	74 68                	je     9fb <_Z7multSEQR6MatrixRKS_S2_ii+0x9db>
 993:	0f b6 45 43          	movzbl 0x43(%rbp),%eax
 997:	0f be f0             	movsbl %al,%esi
 99a:	bf 00 00 00 00       	mov    $0x0,%edi
 99f:	e8 00 00 00 00       	callq  9a4 <_Z7multSEQR6MatrixRKS_S2_ii+0x984>
 9a4:	48 89 c7             	mov    %rax,%rdi
 9a7:	e8 00 00 00 00       	callq  9ac <_Z7multSEQR6MatrixRKS_S2_ii+0x98c>
 9ac:	48 8b 84 24 d8 00 00 	mov    0xd8(%rsp),%rax
 9b3:	00 
 9b4:	64 48 33 04 25 28 00 	xor    %fs:0x28,%rax
 9bb:	00 00 
 9bd:	0f 85 0f 06 00 00    	jne    fd2 <_Z7multSEQR6MatrixRKS_S2_ii+0xfb2>
 9c3:	48 81 c4 e8 00 00 00 	add    $0xe8,%rsp
 9ca:	5b                   	pop    %rbx
 9cb:	5d                   	pop    %rbp
 9cc:	41 5c                	pop    %r12
 9ce:	41 5d                	pop    %r13
 9d0:	41 5e                	pop    %r14
 9d2:	41 5f                	pop    %r15
 9d4:	c3                   	retq   
 9d5:	0f 1f 00             	nopl   (%rax)
 9d8:	0f 28 c2             	movaps %xmm2,%xmm0
 9db:	e9 b5 f9 ff ff       	jmpq   395 <_Z7multSEQR6MatrixRKS_S2_ii+0x375>
 9e0:	48 89 df             	mov    %rbx,%rdi
 9e3:	e8 00 00 00 00       	callq  9e8 <_Z7multSEQR6MatrixRKS_S2_ii+0x9c8>
 9e8:	48 8b 03             	mov    (%rbx),%rax
 9eb:	be 0a 00 00 00       	mov    $0xa,%esi
 9f0:	48 89 df             	mov    %rbx,%rdi
 9f3:	ff 50 30             	callq  *0x30(%rax)
 9f6:	e9 f4 f6 ff ff       	jmpq   ef <_Z7multSEQR6MatrixRKS_S2_ii+0xcf>
 9fb:	48 89 ef             	mov    %rbp,%rdi
 9fe:	e8 00 00 00 00       	callq  a03 <_Z7multSEQR6MatrixRKS_S2_ii+0x9e3>
 a03:	48 8b 45 00          	mov    0x0(%rbp),%rax
 a07:	be 0a 00 00 00       	mov    $0xa,%esi
 a0c:	48 89 ef             	mov    %rbp,%rdi
 a0f:	ff 50 30             	callq  *0x30(%rax)
 a12:	e9 80 ff ff ff       	jmpq   997 <_Z7multSEQR6MatrixRKS_S2_ii+0x977>
 a17:	4c 89 f7             	mov    %r14,%rdi
 a1a:	e8 00 00 00 00       	callq  a1f <_Z7multSEQR6MatrixRKS_S2_ii+0x9ff>
 a1f:	4d 8b 2e             	mov    (%r14),%r13
 a22:	be 0a 00 00 00       	mov    $0xa,%esi
 a27:	4c 89 f7             	mov    %r14,%rdi
 a2a:	41 ff 55 30          	callq  *0x30(%r13)
 a2e:	e9 18 ff ff ff       	jmpq   94b <_Z7multSEQR6MatrixRKS_S2_ii+0x92b>
 a33:	48 89 ef             	mov    %rbp,%rdi
 a36:	e8 00 00 00 00       	callq  a3b <_Z7multSEQR6MatrixRKS_S2_ii+0xa1b>
 a3b:	48 8b 45 00          	mov    0x0(%rbp),%rax
 a3f:	be 0a 00 00 00       	mov    $0xa,%esi
 a44:	48 89 ef             	mov    %rbp,%rdi
 a47:	ff 50 30             	callq  *0x30(%rax)
 a4a:	e9 3b fe ff ff       	jmpq   88a <_Z7multSEQR6MatrixRKS_S2_ii+0x86a>
 a4f:	ba 12 00 00 00       	mov    $0x12,%edx
 a54:	be 00 00 00 00       	mov    $0x0,%esi
 a59:	bf 00 00 00 00       	mov    $0x0,%edi
 a5e:	e8 00 00 00 00       	callq  a63 <_Z7multSEQR6MatrixRKS_S2_ii+0xa43>
 a63:	f2 0f 10 44 24 40    	movsd  0x40(%rsp),%xmm0
 a69:	4c 8b 1d 00 00 00 00 	mov    0x0(%rip),%r11        # a70 <_Z7multSEQR6MatrixRKS_S2_ii+0xa50>
 a70:	bf 00 00 00 00       	mov    $0x0,%edi
 a75:	f2 0f 5e 44 24 38    	divsd  0x38(%rsp),%xmm0
 a7b:	49 8b 73 e8          	mov    -0x18(%r11),%rsi
 a7f:	4c 89 be 00 00 00 00 	mov    %r15,0x0(%rsi)
 a86:	49 8b 6b e8          	mov    -0x18(%r11),%rbp
 a8a:	48 c7 85 00 00 00 00 	movq   $0x4,0x0(%rbp)
 a91:	04 00 00 00 
 a95:	49 8b 5b e8          	mov    -0x18(%r11),%rbx
 a99:	8b 8b 00 00 00 00    	mov    0x0(%rbx),%ecx
 a9f:	81 e1 fb fe ff ff    	and    $0xfffffefb,%ecx
 aa5:	83 c9 04             	or     $0x4,%ecx
 aa8:	89 8b 00 00 00 00    	mov    %ecx,0x0(%rbx)
 aae:	e8 00 00 00 00       	callq  ab3 <_Z7multSEQR6MatrixRKS_S2_ii+0xa93>
 ab3:	4c 8b 08             	mov    (%rax),%r9
 ab6:	49 89 c4             	mov    %rax,%r12
 ab9:	4d 8b 51 e8          	mov    -0x18(%r9),%r10
 abd:	4e 8b b4 10 f0 00 00 	mov    0xf0(%rax,%r10,1),%r14
 ac4:	00 
 ac5:	4d 85 f6             	test   %r14,%r14
 ac8:	0f 84 09 05 00 00    	je     fd7 <_Z7multSEQR6MatrixRKS_S2_ii+0xfb7>
 ace:	41 80 7e 38 00       	cmpb   $0x0,0x38(%r14)
 ad3:	0f 84 dd 04 00 00    	je     fb6 <_Z7multSEQR6MatrixRKS_S2_ii+0xf96>
 ad9:	41 0f b6 46 43       	movzbl 0x43(%r14),%eax
 ade:	4c 89 e7             	mov    %r12,%rdi
 ae1:	0f be f0             	movsbl %al,%esi
 ae4:	e8 00 00 00 00       	callq  ae9 <_Z7multSEQR6MatrixRKS_S2_ii+0xac9>
 ae9:	48 89 c7             	mov    %rax,%rdi
 aec:	e8 00 00 00 00       	callq  af1 <_Z7multSEQR6MatrixRKS_S2_ii+0xad1>
 af1:	e9 57 fd ff ff       	jmpq   84d <_Z7multSEQR6MatrixRKS_S2_ii+0x82d>
 af6:	4c 89 ef             	mov    %r13,%rdi
 af9:	e8 00 00 00 00       	callq  afe <_Z7multSEQR6MatrixRKS_S2_ii+0xade>
 afe:	49 8b 45 00          	mov    0x0(%r13),%rax
 b02:	be 0a 00 00 00       	mov    $0xa,%esi
 b07:	4c 89 ef             	mov    %r13,%rdi
 b0a:	ff 50 30             	callq  *0x30(%rax)
 b0d:	e9 12 fd ff ff       	jmpq   824 <_Z7multSEQR6MatrixRKS_S2_ii+0x804>
 b12:	48 89 ef             	mov    %rbp,%rdi
 b15:	e8 00 00 00 00       	callq  b1a <_Z7multSEQR6MatrixRKS_S2_ii+0xafa>
 b1a:	48 8b 5d 00          	mov    0x0(%rbp),%rbx
 b1e:	be 0a 00 00 00       	mov    $0xa,%esi
 b23:	48 89 ef             	mov    %rbp,%rdi
 b26:	ff 53 30             	callq  *0x30(%rbx)
 b29:	e9 43 fc ff ff       	jmpq   771 <_Z7multSEQR6MatrixRKS_S2_ii+0x751>
 b2e:	48 89 df             	mov    %rbx,%rdi
 b31:	e8 00 00 00 00       	callq  b36 <_Z7multSEQR6MatrixRKS_S2_ii+0xb16>
 b36:	4c 8b 23             	mov    (%rbx),%r12
 b39:	be 0a 00 00 00       	mov    $0xa,%esi
 b3e:	48 89 df             	mov    %rbx,%rdi
 b41:	41 ff 54 24 30       	callq  *0x30(%r12)
 b46:	e9 6f fb ff ff       	jmpq   6ba <_Z7multSEQR6MatrixRKS_S2_ii+0x69a>
 b4b:	4c 89 ef             	mov    %r13,%rdi
 b4e:	e8 00 00 00 00       	callq  b53 <_Z7multSEQR6MatrixRKS_S2_ii+0xb33>
 b53:	49 8b 6d 00          	mov    0x0(%r13),%rbp
 b57:	be 0a 00 00 00       	mov    $0xa,%esi
 b5c:	4c 89 ef             	mov    %r13,%rdi
 b5f:	ff 55 30             	callq  *0x30(%rbp)
 b62:	e9 03 fb ff ff       	jmpq   66a <_Z7multSEQR6MatrixRKS_S2_ii+0x64a>
 b67:	4c 89 e7             	mov    %r12,%rdi
 b6a:	e8 00 00 00 00       	callq  b6f <_Z7multSEQR6MatrixRKS_S2_ii+0xb4f>
 b6f:	49 8b 04 24          	mov    (%r12),%rax
 b73:	be 0a 00 00 00       	mov    $0xa,%esi
 b78:	4c 89 e7             	mov    %r12,%rdi
 b7b:	ff 50 30             	callq  *0x30(%rax)
 b7e:	e9 87 fa ff ff       	jmpq   60a <_Z7multSEQR6MatrixRKS_S2_ii+0x5ea>
 b83:	4c 89 ff             	mov    %r15,%rdi
 b86:	e8 00 00 00 00       	callq  b8b <_Z7multSEQR6MatrixRKS_S2_ii+0xb6b>
 b8b:	49 8b 07             	mov    (%r15),%rax
 b8e:	be 0a 00 00 00       	mov    $0xa,%esi
 b93:	4c 89 ff             	mov    %r15,%rdi
 b96:	ff 50 30             	callq  *0x30(%rax)
 b99:	e9 d2 f8 ff ff       	jmpq   470 <_Z7multSEQR6MatrixRKS_S2_ii+0x450>
 b9e:	4c 89 f7             	mov    %r14,%rdi
 ba1:	e8 00 00 00 00       	callq  ba6 <_Z7multSEQR6MatrixRKS_S2_ii+0xb86>
 ba6:	4d 8b 2e             	mov    (%r14),%r13
 ba9:	be 0a 00 00 00       	mov    $0xa,%esi
 bae:	4c 89 f7             	mov    %r14,%rdi
 bb1:	41 ff 55 30          	callq  *0x30(%r13)
 bb5:	e9 62 f8 ff ff       	jmpq   41c <_Z7multSEQR6MatrixRKS_S2_ii+0x3fc>
 bba:	48 89 ef             	mov    %rbp,%rdi
 bbd:	e8 00 00 00 00       	callq  bc2 <_Z7multSEQR6MatrixRKS_S2_ii+0xba2>
 bc2:	4c 8b 65 00          	mov    0x0(%rbp),%r12
 bc6:	be 0a 00 00 00       	mov    $0xa,%esi
 bcb:	48 89 ef             	mov    %rbp,%rdi
 bce:	41 ff 54 24 30       	callq  *0x30(%r12)
 bd3:	e9 0a f9 ff ff       	jmpq   4e2 <_Z7multSEQR6MatrixRKS_S2_ii+0x4c2>
 bd8:	44 8b 5c 24 60       	mov    0x60(%rsp),%r11d
 bdd:	45 85 db             	test   %r11d,%r11d
 be0:	0f 84 df f7 ff ff    	je     3c5 <_Z7multSEQR6MatrixRKS_S2_ii+0x3a5>
 be6:	44 89 7c 24 6c       	mov    %r15d,0x6c(%rsp)
 beb:	44 8b 7c 24 6c       	mov    0x6c(%rsp),%r15d
 bf0:	31 db                	xor    %ebx,%ebx
 bf2:	31 ed                	xor    %ebp,%ebp
 bf4:	c7 44 24 64 00 00 00 	movl   $0x0,0x64(%rsp)
 bfb:	00 
 bfc:	c7 44 24 68 00 00 00 	movl   $0x0,0x68(%rsp)
 c03:	00 
 c04:	45 85 ff             	test   %r15d,%r15d
 c07:	0f 84 22 03 00 00    	je     f2f <_Z7multSEQR6MatrixRKS_S2_ii+0xf0f>
 c0d:	0f 1f 00             	nopl   (%rax)
 c10:	44 8b 54 24 6c       	mov    0x6c(%rsp),%r10d
 c15:	44 8b 7c 24 64       	mov    0x64(%rsp),%r15d
 c1a:	46 8d 24 2b          	lea    (%rbx,%r13,1),%r12d
 c1e:	44 89 6c 24 38       	mov    %r13d,0x38(%rsp)
 c23:	41 01 ea             	add    %ebp,%r10d
 c26:	44 89 54 24 48       	mov    %r10d,0x48(%rsp)
 c2b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)
 c30:	44 8b 6c 24 38       	mov    0x38(%rsp),%r13d
 c35:	0f 57 d2             	xorps  %xmm2,%xmm2
 c38:	45 85 ed             	test   %r13d,%r13d
 c3b:	0f 84 15 02 00 00    	je     e56 <_Z7multSEQR6MatrixRKS_S2_ii+0xe36>
 c41:	48 8b 7c 24 58       	mov    0x58(%rsp),%rdi
 c46:	4c 8b 44 24 50       	mov    0x50(%rsp),%r8
 c4b:	89 de                	mov    %ebx,%esi
 c4d:	42 8d 04 3b          	lea    (%rbx,%r15,1),%eax
 c51:	89 d9                	mov    %ebx,%ecx
 c53:	44 8d 5b 01          	lea    0x1(%rbx),%r11d
 c57:	f7 d1                	not    %ecx
 c59:	0f 57 c0             	xorps  %xmm0,%xmm0
 c5c:	4c 8b 77 10          	mov    0x10(%rdi),%r14
 c60:	49 8b 50 10          	mov    0x10(%r8),%rdx
 c64:	44 01 e1             	add    %r12d,%ecx
 c67:	83 e1 07             	and    $0x7,%ecx
 c6a:	45 39 e3             	cmp    %r12d,%r11d
 c6d:	f3 41 0f 10 14 b6    	movss  (%r14,%rsi,4),%xmm2
 c73:	f3 0f 59 14 82       	mulss  (%rdx,%rax,4),%xmm2
 c78:	f3 0f 58 d0          	addss  %xmm0,%xmm2
 c7c:	0f 84 d4 01 00 00    	je     e56 <_Z7multSEQR6MatrixRKS_S2_ii+0xe36>
 c82:	85 c9                	test   %ecx,%ecx
 c84:	0f 84 ed 00 00 00    	je     d77 <_Z7multSEQR6MatrixRKS_S2_ii+0xd57>
 c8a:	83 f9 01             	cmp    $0x1,%ecx
 c8d:	0f 84 bf 00 00 00    	je     d52 <_Z7multSEQR6MatrixRKS_S2_ii+0xd32>
 c93:	83 f9 02             	cmp    $0x2,%ecx
 c96:	0f 84 9a 00 00 00    	je     d36 <_Z7multSEQR6MatrixRKS_S2_ii+0xd16>
 c9c:	83 f9 03             	cmp    $0x3,%ecx
 c9f:	90                   	nop
 ca0:	74 79                	je     d1b <_Z7multSEQR6MatrixRKS_S2_ii+0xcfb>
 ca2:	83 f9 04             	cmp    $0x4,%ecx
 ca5:	74 5a                	je     d01 <_Z7multSEQR6MatrixRKS_S2_ii+0xce1>
 ca7:	83 f9 05             	cmp    $0x5,%ecx
 caa:	74 3b                	je     ce7 <_Z7multSEQR6MatrixRKS_S2_ii+0xcc7>
 cac:	83 f9 06             	cmp    $0x6,%ecx
 caf:	90                   	nop
 cb0:	74 1a                	je     ccc <_Z7multSEQR6MatrixRKS_S2_ii+0xcac>
 cb2:	45 89 d9             	mov    %r11d,%r9d
 cb5:	45 01 fb             	add    %r15d,%r11d
 cb8:	f3 43 0f 10 1c 8e    	movss  (%r14,%r9,4),%xmm3
 cbe:	f3 42 0f 59 1c 9a    	mulss  (%rdx,%r11,4),%xmm3
 cc4:	44 8d 5b 02          	lea    0x2(%rbx),%r11d
 cc8:	f3 0f 58 d3          	addss  %xmm3,%xmm2
 ccc:	45 89 da             	mov    %r11d,%r10d
 ccf:	47 8d 2c 3b          	lea    (%r11,%r15,1),%r13d
 cd3:	41 83 c3 01          	add    $0x1,%r11d
 cd7:	f3 43 0f 10 24 96    	movss  (%r14,%r10,4),%xmm4
 cdd:	f3 42 0f 59 24 aa    	mulss  (%rdx,%r13,4),%xmm4
 ce3:	f3 0f 58 d4          	addss  %xmm4,%xmm2
 ce7:	45 89 d8             	mov    %r11d,%r8d
 cea:	43 8d 3c 3b          	lea    (%r11,%r15,1),%edi
 cee:	41 83 c3 01          	add    $0x1,%r11d
 cf2:	f3 43 0f 10 2c 86    	movss  (%r14,%r8,4),%xmm5
 cf8:	f3 0f 59 2c ba       	mulss  (%rdx,%rdi,4),%xmm5
 cfd:	f3 0f 58 d5          	addss  %xmm5,%xmm2
 d01:	44 89 d9             	mov    %r11d,%ecx
 d04:	43 8d 34 3b          	lea    (%r11,%r15,1),%esi
 d08:	41 83 c3 01          	add    $0x1,%r11d
 d0c:	f3 41 0f 10 34 8e    	movss  (%r14,%rcx,4),%xmm6
 d12:	f3 0f 59 34 b2       	mulss  (%rdx,%rsi,4),%xmm6
 d17:	f3 0f 58 d6          	addss  %xmm6,%xmm2
 d1b:	44 89 d8             	mov    %r11d,%eax
 d1e:	47 8d 0c 3b          	lea    (%r11,%r15,1),%r9d
 d22:	41 83 c3 01          	add    $0x1,%r11d
 d26:	f3 41 0f 10 3c 86    	movss  (%r14,%rax,4),%xmm7
 d2c:	f3 42 0f 59 3c 8a    	mulss  (%rdx,%r9,4),%xmm7
 d32:	f3 0f 58 d7          	addss  %xmm7,%xmm2
 d36:	45 89 da             	mov    %r11d,%r10d
 d39:	47 8d 2c 3b          	lea    (%r11,%r15,1),%r13d
 d3d:	41 83 c3 01          	add    $0x1,%r11d
 d41:	f3 47 0f 10 04 96    	movss  (%r14,%r10,4),%xmm8
 d47:	f3 46 0f 59 04 aa    	mulss  (%rdx,%r13,4),%xmm8
 d4d:	f3 41 0f 58 d0       	addss  %xmm8,%xmm2
 d52:	45 89 d8             	mov    %r11d,%r8d
 d55:	43 8d 3c 3b          	lea    (%r11,%r15,1),%edi
 d59:	41 83 c3 01          	add    $0x1,%r11d
 d5d:	f3 47 0f 10 0c 86    	movss  (%r14,%r8,4),%xmm9
 d63:	45 39 e3             	cmp    %r12d,%r11d
 d66:	f3 44 0f 59 0c ba    	mulss  (%rdx,%rdi,4),%xmm9
 d6c:	f3 41 0f 58 d1       	addss  %xmm9,%xmm2
 d71:	0f 84 df 00 00 00    	je     e56 <_Z7multSEQR6MatrixRKS_S2_ii+0xe36>
 d77:	44 89 de             	mov    %r11d,%esi
 d7a:	43 8d 0c 3b          	lea    (%r11,%r15,1),%ecx
 d7e:	41 8d 43 01          	lea    0x1(%r11),%eax
 d82:	f3 45 0f 10 14 b6    	movss  (%r14,%rsi,4),%xmm10
 d88:	45 8d 53 02          	lea    0x2(%r11),%r10d
 d8c:	f3 44 0f 59 14 8a    	mulss  (%rdx,%rcx,4),%xmm10
 d92:	41 89 c1             	mov    %eax,%r9d
 d95:	44 01 f8             	add    %r15d,%eax
 d98:	f3 47 0f 10 1c 8e    	movss  (%r14,%r9,4),%xmm11
 d9e:	45 89 d5             	mov    %r10d,%r13d
 da1:	f3 44 0f 59 1c 82    	mulss  (%rdx,%rax,4),%xmm11
 da7:	45 01 fa             	add    %r15d,%r10d
 daa:	45 8d 43 03          	lea    0x3(%r11),%r8d
 dae:	f3 47 0f 10 24 ae    	movss  (%r14,%r13,4),%xmm12
 db4:	41 8d 73 04          	lea    0x4(%r11),%esi
 db8:	f3 46 0f 59 24 92    	mulss  (%rdx,%r10,4),%xmm12
 dbe:	44 89 c7             	mov    %r8d,%edi
 dc1:	45 01 f8             	add    %r15d,%r8d
 dc4:	f3 45 0f 10 2c be    	movss  (%r14,%rdi,4),%xmm13
 dca:	89 f1                	mov    %esi,%ecx
 dcc:	f3 41 0f 58 d2       	addss  %xmm10,%xmm2
 dd1:	f3 46 0f 59 2c 82    	mulss  (%rdx,%r8,4),%xmm13
 dd7:	44 01 fe             	add    %r15d,%esi
 dda:	f3 45 0f 10 34 8e    	movss  (%r14,%rcx,4),%xmm14
 de0:	41 8d 43 05          	lea    0x5(%r11),%eax
 de4:	f3 44 0f 59 34 b2    	mulss  (%rdx,%rsi,4),%xmm14
 dea:	45 8d 53 06          	lea    0x6(%r11),%r10d
 dee:	45 8d 43 07          	lea    0x7(%r11),%r8d
 df2:	41 89 c1             	mov    %eax,%r9d
 df5:	44 01 f8             	add    %r15d,%eax
 df8:	41 83 c3 08          	add    $0x8,%r11d
 dfc:	f3 41 0f 58 d3       	addss  %xmm11,%xmm2
 e01:	f3 47 0f 10 3c 8e    	movss  (%r14,%r9,4),%xmm15
 e07:	f3 44 0f 59 3c 82    	mulss  (%rdx,%rax,4),%xmm15
 e0d:	45 89 d5             	mov    %r10d,%r13d
 e10:	45 01 fa             	add    %r15d,%r10d
 e13:	f3 43 0f 10 0c ae    	movss  (%r14,%r13,4),%xmm1
 e19:	44 89 c7             	mov    %r8d,%edi
 e1c:	f3 42 0f 59 0c 92    	mulss  (%rdx,%r10,4),%xmm1
 e22:	45 01 f8             	add    %r15d,%r8d
 e25:	45 39 e3             	cmp    %r12d,%r11d
 e28:	f3 41 0f 58 d4       	addss  %xmm12,%xmm2
 e2d:	f3 41 0f 10 04 be    	movss  (%r14,%rdi,4),%xmm0
 e33:	f3 42 0f 59 04 82    	mulss  (%rdx,%r8,4),%xmm0
 e39:	f3 41 0f 58 d5       	addss  %xmm13,%xmm2
 e3e:	f3 41 0f 58 d6       	addss  %xmm14,%xmm2
 e43:	f3 41 0f 58 d7       	addss  %xmm15,%xmm2
 e48:	f3 0f 58 d1          	addss  %xmm1,%xmm2
 e4c:	f3 0f 58 d0          	addss  %xmm0,%xmm2
 e50:	0f 85 21 ff ff ff    	jne    d77 <_Z7multSEQR6MatrixRKS_S2_ii+0xd57>
 e56:	41 89 ee             	mov    %ebp,%r14d
 e59:	bf 00 00 00 00       	mov    $0x0,%edi
 e5e:	f3 0f 11 14 24       	movss  %xmm2,(%rsp)
 e63:	4c 89 f6             	mov    %r14,%rsi
 e66:	e8 00 00 00 00       	callq  e6b <_Z7multSEQR6MatrixRKS_S2_ii+0xe4b>
 e6b:	ba 02 00 00 00       	mov    $0x2,%edx
 e70:	be 00 00 00 00       	mov    $0x0,%esi
 e75:	48 89 c7             	mov    %rax,%rdi
 e78:	49 89 c5             	mov    %rax,%r13
 e7b:	e8 00 00 00 00       	callq  e80 <_Z7multSEQR6MatrixRKS_S2_ii+0xe60>
 e80:	f3 0f 10 14 24       	movss  (%rsp),%xmm2
 e85:	4c 89 ef             	mov    %r13,%rdi
 e88:	0f 14 d2             	unpcklps %xmm2,%xmm2
 e8b:	0f 5a c2             	cvtps2pd %xmm2,%xmm0
 e8e:	e8 00 00 00 00       	callq  e93 <_Z7multSEQR6MatrixRKS_S2_ii+0xe73>
 e93:	48 8b 10             	mov    (%rax),%rdx
 e96:	49 89 c1             	mov    %rax,%r9
 e99:	f3 0f 10 1c 24       	movss  (%rsp),%xmm3
 e9e:	4c 8b 5a e8          	mov    -0x18(%rdx),%r11
 ea2:	4a 8b b4 18 f0 00 00 	mov    0xf0(%rax,%r11,1),%rsi
 ea9:	00 
 eaa:	48 85 f6             	test   %rsi,%rsi
 ead:	0f 84 24 01 00 00    	je     fd7 <_Z7multSEQR6MatrixRKS_S2_ii+0xfb7>
 eb3:	80 7e 38 00          	cmpb   $0x0,0x38(%rsi)
 eb7:	74 7f                	je     f38 <_Z7multSEQR6MatrixRKS_S2_ii+0xf18>
 eb9:	0f b6 46 43          	movzbl 0x43(%rsi),%eax
 ebd:	4c 89 cf             	mov    %r9,%rdi
 ec0:	0f be f0             	movsbl %al,%esi
 ec3:	f3 0f 11 1c 24       	movss  %xmm3,(%rsp)
 ec8:	e8 00 00 00 00       	callq  ecd <_Z7multSEQR6MatrixRKS_S2_ii+0xead>
 ecd:	48 89 c7             	mov    %rax,%rdi
 ed0:	83 c5 01             	add    $0x1,%ebp
 ed3:	e8 00 00 00 00       	callq  ed8 <_Z7multSEQR6MatrixRKS_S2_ii+0xeb8>
 ed8:	4c 8b 54 24 40       	mov    0x40(%rsp),%r10
 edd:	44 03 7c 24 38       	add    0x38(%rsp),%r15d
 ee2:	3b 6c 24 48          	cmp    0x48(%rsp),%ebp
 ee6:	f3 0f 10 24 24       	movss  (%rsp),%xmm4
 eeb:	4d 8b 42 10          	mov    0x10(%r10),%r8
 eef:	f3 43 0f 11 24 b0    	movss  %xmm4,(%r8,%r14,4)
 ef5:	0f 85 35 fd ff ff    	jne    c30 <_Z7multSEQR6MatrixRKS_S2_ii+0xc10>
 efb:	44 8b 6c 24 38       	mov    0x38(%rsp),%r13d
 f00:	83 44 24 68 01       	addl   $0x1,0x68(%rsp)
 f05:	44 29 6c 24 64       	sub    %r13d,0x64(%rsp)
 f0a:	44 01 eb             	add    %r13d,%ebx
 f0d:	44 8b 7c 24 60       	mov    0x60(%rsp),%r15d
 f12:	44 39 7c 24 68       	cmp    %r15d,0x68(%rsp)
 f17:	8b 6c 24 48          	mov    0x48(%rsp),%ebp
 f1b:	0f 84 a4 f4 ff ff    	je     3c5 <_Z7multSEQR6MatrixRKS_S2_ii+0x3a5>
 f21:	44 8b 7c 24 6c       	mov    0x6c(%rsp),%r15d
 f26:	45 85 ff             	test   %r15d,%r15d
 f29:	0f 85 e1 fc ff ff    	jne    c10 <_Z7multSEQR6MatrixRKS_S2_ii+0xbf0>
 f2f:	89 6c 24 48          	mov    %ebp,0x48(%rsp)
 f33:	eb cb                	jmp    f00 <_Z7multSEQR6MatrixRKS_S2_ii+0xee0>
 f35:	0f 1f 00             	nopl   (%rax)
 f38:	48 89 f7             	mov    %rsi,%rdi
 f3b:	f3 0f 11 1c 24       	movss  %xmm3,(%rsp)
 f40:	48 89 74 24 20       	mov    %rsi,0x20(%rsp)
 f45:	48 89 44 24 18       	mov    %rax,0x18(%rsp)
 f4a:	e8 00 00 00 00       	callq  f4f <_Z7multSEQR6MatrixRKS_S2_ii+0xf2f>
 f4f:	48 8b 4c 24 20       	mov    0x20(%rsp),%rcx
 f54:	be 0a 00 00 00       	mov    $0xa,%esi
 f59:	48 8b 01             	mov    (%rcx),%rax
 f5c:	48 89 cf             	mov    %rcx,%rdi
 f5f:	ff 50 30             	callq  *0x30(%rax)
 f62:	4c 8b 4c 24 18       	mov    0x18(%rsp),%r9
 f67:	f3 0f 10 1c 24       	movss  (%rsp),%xmm3
 f6c:	e9 4c ff ff ff       	jmpq   ebd <_Z7multSEQR6MatrixRKS_S2_ii+0xe9d>
 f71:	49 d1 e9             	shr    %r9
 f74:	f2 49 0f 2a f1       	cvtsi2sd %r9,%xmm6
 f79:	f2 0f 58 f6          	addsd  %xmm6,%xmm6
 f7d:	f2 0f 11 74 24 48    	movsd  %xmm6,0x48(%rsp)
 f83:	e9 a2 f5 ff ff       	jmpq   52a <_Z7multSEQR6MatrixRKS_S2_ii+0x50a>
 f88:	48 89 ee             	mov    %rbp,%rsi
 f8b:	4c 89 c7             	mov    %r8,%rdi
 f8e:	4c 29 d6             	sub    %r10,%rsi
 f91:	48 83 c7 10          	add    $0x10,%rdi
 f95:	e8 00 00 00 00       	callq  f9a <_Z7multSEQR6MatrixRKS_S2_ii+0xf7a>
 f9a:	e9 13 f1 ff ff       	jmpq   b2 <_Z7multSEQR6MatrixRKS_S2_ii+0x92>
 f9f:	48 8b 7c 24 50       	mov    0x50(%rsp),%rdi
 fa4:	8b 36                	mov    (%rsi),%esi
 fa6:	44 8b 3f             	mov    (%rdi),%r15d
 fa9:	89 74 24 60          	mov    %esi,0x60(%rsp)
 fad:	44 8b 6f 04          	mov    0x4(%rdi),%r13d
 fb1:	e9 c7 f0 ff ff       	jmpq   7d <_Z7multSEQR6MatrixRKS_S2_ii+0x5d>
 fb6:	4c 89 f7             	mov    %r14,%rdi
 fb9:	e8 00 00 00 00       	callq  fbe <_Z7multSEQR6MatrixRKS_S2_ii+0xf9e>
 fbe:	4d 8b 2e             	mov    (%r14),%r13
 fc1:	be 0a 00 00 00       	mov    $0xa,%esi
 fc6:	4c 89 f7             	mov    %r14,%rdi
 fc9:	41 ff 55 30          	callq  *0x30(%r13)
 fcd:	e9 0c fb ff ff       	jmpq   ade <_Z7multSEQR6MatrixRKS_S2_ii+0xabe>
 fd2:	e8 00 00 00 00       	callq  fd7 <_Z7multSEQR6MatrixRKS_S2_ii+0xfb7>
 fd7:	e8 00 00 00 00       	callq  fdc <_Z7multSEQR6MatrixRKS_S2_ii+0xfbc>

Disassembly of section .text._ZNSt6vectorIfSaIfEE17_M_default_appendEm:

0000000000000000 <_ZNSt6vectorIfSaIfEE17_M_default_appendEm>:
   0:	48 89 5c 24 d8       	mov    %rbx,-0x28(%rsp)
   5:	48 89 6c 24 e0       	mov    %rbp,-0x20(%rsp)
   a:	48 89 f3             	mov    %rsi,%rbx
   d:	4c 89 64 24 e8       	mov    %r12,-0x18(%rsp)
  12:	4c 89 6c 24 f0       	mov    %r13,-0x10(%rsp)
  17:	4c 89 74 24 f8       	mov    %r14,-0x8(%rsp)
  1c:	48 83 ec 28          	sub    $0x28,%rsp
  20:	48 85 f6             	test   %rsi,%rsi
  23:	0f 84 14 01 00 00    	je     13d <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0x13d>
  29:	4c 8b 5f 08          	mov    0x8(%rdi),%r11
  2d:	48 8b 47 10          	mov    0x10(%rdi),%rax
  31:	48 89 fd             	mov    %rdi,%rbp
  34:	4c 29 d8             	sub    %r11,%rax
  37:	48 c1 f8 02          	sar    $0x2,%rax
  3b:	48 39 c6             	cmp    %rax,%rsi
  3e:	0f 87 1c 01 00 00    	ja     160 <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0x160>
  44:	48 8d 7e ff          	lea    -0x1(%rsi),%rdi
  48:	49 8d 73 04          	lea    0x4(%r11),%rsi
  4c:	41 c7 03 00 00 00 00 	movl   $0x0,(%r11)
  53:	49 89 f8             	mov    %rdi,%r8
  56:	41 83 e0 07          	and    $0x7,%r8d
  5a:	48 85 ff             	test   %rdi,%rdi
  5d:	0f 84 d2 00 00 00    	je     135 <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0x135>
  63:	4d 85 c0             	test   %r8,%r8
  66:	0f 84 88 00 00 00    	je     f4 <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0xf4>
  6c:	49 83 f8 01          	cmp    $0x1,%r8
  70:	74 72                	je     e4 <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0xe4>
  72:	49 83 f8 02          	cmp    $0x2,%r8
  76:	74 5e                	je     d6 <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0xd6>
  78:	49 83 f8 03          	cmp    $0x3,%r8
  7c:	74 4a                	je     c8 <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0xc8>
  7e:	49 83 f8 04          	cmp    $0x4,%r8
  82:	74 36                	je     ba <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0xba>
  84:	49 83 f8 05          	cmp    $0x5,%r8
  88:	74 22                	je     ac <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0xac>
  8a:	49 83 f8 06          	cmp    $0x6,%r8
  8e:	74 0e                	je     9e <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0x9e>
  90:	c7 06 00 00 00 00    	movl   $0x0,(%rsi)
  96:	48 83 ef 01          	sub    $0x1,%rdi
  9a:	48 83 c6 04          	add    $0x4,%rsi
  9e:	c7 06 00 00 00 00    	movl   $0x0,(%rsi)
  a4:	48 83 ef 01          	sub    $0x1,%rdi
  a8:	48 83 c6 04          	add    $0x4,%rsi
  ac:	c7 06 00 00 00 00    	movl   $0x0,(%rsi)
  b2:	48 83 ef 01          	sub    $0x1,%rdi
  b6:	48 83 c6 04          	add    $0x4,%rsi
  ba:	c7 06 00 00 00 00    	movl   $0x0,(%rsi)
  c0:	48 83 ef 01          	sub    $0x1,%rdi
  c4:	48 83 c6 04          	add    $0x4,%rsi
  c8:	c7 06 00 00 00 00    	movl   $0x0,(%rsi)
  ce:	48 83 ef 01          	sub    $0x1,%rdi
  d2:	48 83 c6 04          	add    $0x4,%rsi
  d6:	c7 06 00 00 00 00    	movl   $0x0,(%rsi)
  dc:	48 83 ef 01          	sub    $0x1,%rdi
  e0:	48 83 c6 04          	add    $0x4,%rsi
  e4:	c7 06 00 00 00 00    	movl   $0x0,(%rsi)
  ea:	48 83 c6 04          	add    $0x4,%rsi
  ee:	48 83 ef 01          	sub    $0x1,%rdi
  f2:	74 41                	je     135 <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0x135>
  f4:	c7 06 00 00 00 00    	movl   $0x0,(%rsi)
  fa:	c7 46 04 00 00 00 00 	movl   $0x0,0x4(%rsi)
 101:	c7 46 08 00 00 00 00 	movl   $0x0,0x8(%rsi)
 108:	c7 46 0c 00 00 00 00 	movl   $0x0,0xc(%rsi)
 10f:	c7 46 10 00 00 00 00 	movl   $0x0,0x10(%rsi)
 116:	c7 46 14 00 00 00 00 	movl   $0x0,0x14(%rsi)
 11d:	c7 46 18 00 00 00 00 	movl   $0x0,0x18(%rsi)
 124:	c7 46 1c 00 00 00 00 	movl   $0x0,0x1c(%rsi)
 12b:	48 83 c6 20          	add    $0x20,%rsi
 12f:	48 83 ef 08          	sub    $0x8,%rdi
 133:	75 bf                	jne    f4 <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0xf4>
 135:	4d 8d 0c 9b          	lea    (%r11,%rbx,4),%r9
 139:	4c 89 4d 08          	mov    %r9,0x8(%rbp)
 13d:	48 8b 1c 24          	mov    (%rsp),%rbx
 141:	48 8b 6c 24 08       	mov    0x8(%rsp),%rbp
 146:	4c 8b 64 24 10       	mov    0x10(%rsp),%r12
 14b:	4c 8b 6c 24 18       	mov    0x18(%rsp),%r13
 150:	4c 8b 74 24 20       	mov    0x20(%rsp),%r14
 155:	48 83 c4 28          	add    $0x28,%rsp
 159:	c3                   	retq   
 15a:	66 0f 1f 44 00 00    	nopw   0x0(%rax,%rax,1)
 160:	48 8b 3f             	mov    (%rdi),%rdi
 163:	48 be ff ff ff ff ff 	movabs $0x3fffffffffffffff,%rsi
 16a:	ff ff 3f 
 16d:	49 89 f1             	mov    %rsi,%r9
 170:	49 29 fb             	sub    %rdi,%r11
 173:	49 c1 fb 02          	sar    $0x2,%r11
 177:	4d 29 d9             	sub    %r11,%r9
 17a:	4c 89 da             	mov    %r11,%rdx
 17d:	4d 89 d8             	mov    %r11,%r8
 180:	4c 39 cb             	cmp    %r9,%rbx
 183:	0f 87 90 01 00 00    	ja     319 <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0x319>
 189:	49 39 db             	cmp    %rbx,%r11
 18c:	49 89 da             	mov    %rbx,%r10
 18f:	49 c7 c4 fc ff ff ff 	mov    $0xfffffffffffffffc,%r12
 196:	4d 0f 43 d3          	cmovae %r11,%r10
 19a:	4d 01 d3             	add    %r10,%r11
 19d:	0f 83 5d 01 00 00    	jae    300 <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0x300>
 1a3:	4c 89 e7             	mov    %r12,%rdi
 1a6:	e8 00 00 00 00       	callq  1ab <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0x1ab>
 1ab:	48 8b 7d 00          	mov    0x0(%rbp),%rdi
 1af:	48 8b 55 08          	mov    0x8(%rbp),%rdx
 1b3:	49 89 c5             	mov    %rax,%r13
 1b6:	48 29 fa             	sub    %rdi,%rdx
 1b9:	48 c1 fa 02          	sar    $0x2,%rdx
 1bd:	49 89 d0             	mov    %rdx,%r8
 1c0:	48 85 d2             	test   %rdx,%rdx
 1c3:	4e 8d 34 85 00 00 00 	lea    0x0(,%r8,4),%r14
 1ca:	00 
 1cb:	74 12                	je     1df <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0x1df>
 1cd:	48 89 fe             	mov    %rdi,%rsi
 1d0:	4c 89 f2             	mov    %r14,%rdx
 1d3:	4c 89 ef             	mov    %r13,%rdi
 1d6:	e8 00 00 00 00       	callq  1db <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0x1db>
 1db:	48 8b 7d 00          	mov    0x0(%rbp),%rdi
 1df:	48 8d 4b ff          	lea    -0x1(%rbx),%rcx
 1e3:	4b 8d 54 35 00       	lea    0x0(%r13,%r14,1),%rdx
 1e8:	49 89 cb             	mov    %rcx,%r11
 1eb:	48 8d 42 04          	lea    0x4(%rdx),%rax
 1ef:	41 83 e3 07          	and    $0x7,%r11d
 1f3:	48 85 c9             	test   %rcx,%rcx
 1f6:	c7 02 00 00 00 00    	movl   $0x0,(%rdx)
 1fc:	0f 84 d2 00 00 00    	je     2d4 <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0x2d4>
 202:	4d 85 db             	test   %r11,%r11
 205:	0f 84 88 00 00 00    	je     293 <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0x293>
 20b:	49 83 fb 01          	cmp    $0x1,%r11
 20f:	74 72                	je     283 <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0x283>
 211:	49 83 fb 02          	cmp    $0x2,%r11
 215:	74 5e                	je     275 <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0x275>
 217:	49 83 fb 03          	cmp    $0x3,%r11
 21b:	74 4a                	je     267 <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0x267>
 21d:	49 83 fb 04          	cmp    $0x4,%r11
 221:	74 36                	je     259 <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0x259>
 223:	49 83 fb 05          	cmp    $0x5,%r11
 227:	74 22                	je     24b <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0x24b>
 229:	49 83 fb 06          	cmp    $0x6,%r11
 22d:	74 0e                	je     23d <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0x23d>
 22f:	c7 00 00 00 00 00    	movl   $0x0,(%rax)
 235:	48 8d 4b fe          	lea    -0x2(%rbx),%rcx
 239:	48 8d 42 08          	lea    0x8(%rdx),%rax
 23d:	c7 00 00 00 00 00    	movl   $0x0,(%rax)
 243:	48 83 e9 01          	sub    $0x1,%rcx
 247:	48 83 c0 04          	add    $0x4,%rax
 24b:	c7 00 00 00 00 00    	movl   $0x0,(%rax)
 251:	48 83 e9 01          	sub    $0x1,%rcx
 255:	48 83 c0 04          	add    $0x4,%rax
 259:	c7 00 00 00 00 00    	movl   $0x0,(%rax)
 25f:	48 83 e9 01          	sub    $0x1,%rcx
 263:	48 83 c0 04          	add    $0x4,%rax
 267:	c7 00 00 00 00 00    	movl   $0x0,(%rax)
 26d:	48 83 e9 01          	sub    $0x1,%rcx
 271:	48 83 c0 04          	add    $0x4,%rax
 275:	c7 00 00 00 00 00    	movl   $0x0,(%rax)
 27b:	48 83 e9 01          	sub    $0x1,%rcx
 27f:	48 83 c0 04          	add    $0x4,%rax
 283:	c7 00 00 00 00 00    	movl   $0x0,(%rax)
 289:	48 83 c0 04          	add    $0x4,%rax
 28d:	48 83 e9 01          	sub    $0x1,%rcx
 291:	74 41                	je     2d4 <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0x2d4>
 293:	c7 00 00 00 00 00    	movl   $0x0,(%rax)
 299:	c7 40 04 00 00 00 00 	movl   $0x0,0x4(%rax)
 2a0:	c7 40 08 00 00 00 00 	movl   $0x0,0x8(%rax)
 2a7:	c7 40 0c 00 00 00 00 	movl   $0x0,0xc(%rax)
 2ae:	c7 40 10 00 00 00 00 	movl   $0x0,0x10(%rax)
 2b5:	c7 40 14 00 00 00 00 	movl   $0x0,0x14(%rax)
 2bc:	c7 40 18 00 00 00 00 	movl   $0x0,0x18(%rax)
 2c3:	c7 40 1c 00 00 00 00 	movl   $0x0,0x1c(%rax)
 2ca:	48 83 c0 20          	add    $0x20,%rax
 2ce:	48 83 e9 08          	sub    $0x8,%rcx
 2d2:	75 bf                	jne    293 <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0x293>
 2d4:	48 85 ff             	test   %rdi,%rdi
 2d7:	4c 8d 34 9a          	lea    (%rdx,%rbx,4),%r14
 2db:	74 05                	je     2e2 <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0x2e2>
 2dd:	e8 00 00 00 00       	callq  2e2 <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0x2e2>
 2e2:	4d 01 ec             	add    %r13,%r12
 2e5:	4c 89 6d 00          	mov    %r13,0x0(%rbp)
 2e9:	4c 89 75 08          	mov    %r14,0x8(%rbp)
 2ed:	4c 89 65 10          	mov    %r12,0x10(%rbp)
 2f1:	e9 47 fe ff ff       	jmpq   13d <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0x13d>
 2f6:	66 2e 0f 1f 84 00 00 	nopw   %cs:0x0(%rax,%rax,1)
 2fd:	00 00 00 
 300:	49 39 f3             	cmp    %rsi,%r11
 303:	0f 87 9a fe ff ff    	ja     1a3 <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0x1a3>
 309:	4d 85 db             	test   %r11,%r11
 30c:	75 15                	jne    323 <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0x323>
 30e:	45 31 e4             	xor    %r12d,%r12d
 311:	45 31 ed             	xor    %r13d,%r13d
 314:	e9 a7 fe ff ff       	jmpq   1c0 <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0x1c0>
 319:	bf 00 00 00 00       	mov    $0x0,%edi
 31e:	e8 00 00 00 00       	callq  323 <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0x323>
 323:	4e 8d 24 9d 00 00 00 	lea    0x0(,%r11,4),%r12
 32a:	00 
 32b:	e9 73 fe ff ff       	jmpq   1a3 <_ZNSt6vectorIfSaIfEE17_M_default_appendEm+0x1a3>

Disassembly of section .text.startup:

0000000000000000 <_GLOBAL__sub_I__Z7multSEQR6MatrixRKS_S2_ii>:
   0:	48 83 ec 08          	sub    $0x8,%rsp
   4:	bf 00 00 00 00       	mov    $0x0,%edi
   9:	e8 00 00 00 00       	callq  e <_GLOBAL__sub_I__Z7multSEQR6MatrixRKS_S2_ii+0xe>
   e:	ba 00 00 00 00       	mov    $0x0,%edx
  13:	be 00 00 00 00       	mov    $0x0,%esi
  18:	bf 00 00 00 00       	mov    $0x0,%edi
  1d:	e8 00 00 00 00       	callq  22 <_GLOBAL__sub_I__Z7multSEQR6MatrixRKS_S2_ii+0x22>
  22:	66 c7 05 00 00 00 00 	movw   $0xffff,0x0(%rip)        # 2b <_GLOBAL__sub_I__Z7multSEQR6MatrixRKS_S2_ii+0x2b>
  29:	ff ff 
  2b:	66 c7 05 00 00 00 00 	movw   $0xffff,0x0(%rip)        # 34 <_GLOBAL__sub_I__Z7multSEQR6MatrixRKS_S2_ii+0x34>
  32:	ff ff 
  34:	66 c7 05 00 00 00 00 	movw   $0xfffe,0x0(%rip)        # 3d <_GLOBAL__sub_I__Z7multSEQR6MatrixRKS_S2_ii+0x3d>
  3b:	fe ff 
  3d:	66 c7 05 00 00 00 00 	movw   $0xfffe,0x0(%rip)        # 46 <_GLOBAL__sub_I__Z7multSEQR6MatrixRKS_S2_ii+0x46>
  44:	fe ff 
  46:	48 83 c4 08          	add    $0x8,%rsp
  4a:	c3                   	retq   
